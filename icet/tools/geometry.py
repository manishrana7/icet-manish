import numpy as np
from ase import Atoms
import spglib

from icet.core.lattice_site import LatticeSite
from icet.core_py.lattice_site import LatticeSite as LatticeSite_py


def get_scaled_positions(positions, cell, wrap=True, pbc=[True, True, True]):
    '''Get positions in reduced (scaled) coordinates.

    If wrap is True, positions outside the unit cell will be wrapped into
    the cell in the directions with periodic boundary conditions
    such that the scaled coordinates are between zero and one.
    '''

    fractional = np.linalg.solve(cell.T, positions.T).T

    if wrap:
        for i, periodic in enumerate(pbc):
            if periodic:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test.
                fractional[:, i] %= 1.0
                fractional[:, i] %= 1.0

    return fractional


def find_lattice_site_from_position_python(structure, position):
    '''
    Get lattice neighbor from position.

    This is the Python version of
    `structure.findLatticeSiteFromPosition(position)`

    It is slower but kept for debugging and if further development is needed.
    '''

    fractional = np.linalg.solve(structure.cell.T, np.array(position).T).T
    unit_cell_offset = [int(round(x)) for x in fractional]

    residual = np.dot(fractional - unit_cell_offset, structure.cell)
    try:
        index = structure.find_index_of_position(residual)
    except Exception:
        msg = ['error did not find index with pos: {}'.format(residual)]
        msg += ['position in structure are:']
        msg += ['\n' + str(structure.positions)]
        raise Exception(' '.join(msg))

    latNbr = LatticeSite(index, unit_cell_offset)
    return latNbr


def add_vacuum_in_non_pbc(atoms):
    '''
    Add vacuum in non-periodic directions.

    Parameters
    ----------
    atoms : ASE Atoms object
        input structure

    Returns
    -------
    ASE Atoms object
        output structure
    '''

    vacuum_axis = []
    for i, pbc in enumerate(atoms.pbc):
        if not pbc:
            vacuum_axis.append(i)

    if len(vacuum_axis) > 0:
        atoms.center(30, axis=vacuum_axis)
        atoms.wrap()

    return atoms


def get_primitive_structure(atoms, no_idealize=True):
    '''
    Determines primitive structure using spglib.

    Parameters
    ----------
    atoms : ASE Atoms object
        input structure

    Returns
    -------
    ASE Atoms object
        output structure
    '''
    atoms_with_vacuum = add_vacuum_in_non_pbc(atoms)

    lattice, scaled_positions, numbers = spglib.standardize_cell(
        atoms_with_vacuum, to_primitive=True, no_idealize=no_idealize)
    scaled_positions = [np.round(pos, 12) for pos in scaled_positions]
    atoms_prim = Atoms(scaled_positions=scaled_positions,
                       numbers=numbers, cell=lattice, pbc=atoms.pbc)
    atoms_prim.wrap()
    # icet_wrap(atoms_prim)
    # print(atoms_prim.positions)
    # exit(1)
    return atoms_prim


def get_fractional_positions_from_neighbor_list(structure, neighbor_list):
    '''
    Returns the fractional positions in structure from the neighbors in the
    neighbor list.
    '''
    neighbor_positions = []
    fractional_positions = []
    lattice_site = LatticeSite(0, [0, 0, 0])
    for i in range(len(neighbor_list)):
        lattice_site.index = i
        position = structure.get_position(lattice_site)
        neighbor_positions.append(position)
        for neighbor in neighbor_list.get_neighbors(i):
            position = structure.get_position(neighbor)
            neighbor_positions.append(position)
    if len(neighbor_positions) > 0:
        fractional_positions = get_scaled_positions(
            np.array(neighbor_positions),
            structure.cell, wrap=False,
            pbc=structure.pbc)
    return fractional_positions


def get_permutation_matrix(input_configuration,
                           reference_structure,
                           tolerance_cell=0.05,
                           ):
    '''
    Computes and returns the permutation
    matrix that takes the reference cell to the input cell,
    i.e. permutation_matrix * reference_cell = input_cell
    '''

    input_cell = input_configuration.cell

    # obtain the (in general non-integer) transformation matrix
    # connecting the input configuration to the reference structure
    # L = L_p.P --> P = L_p^-1.L
    P = np.dot(input_cell, np.linalg.inv(reference_structure.cell))

    # assert that the transformation matrix does not deviate too
    # strongly from the nearest integer matrix
    if np.linalg.norm(P - np.around(P)) / 9 > tolerance_cell:
        s = 'Failed to map configuration to reference'
        s += 'structure (tolerance_cell exceeded).\n'
        s += 'reference:\n {}\n'.format(reference_structure.cell)
        s += 'input:\n {}\n'.format(input_configuration.cell)
        s += 'input_cell:\n {}\n'.format(input_cell)
        s += 'P:\n {}\n'.format(P)
        s += 'P_round:\n {}\n'.format(np.around(P))
        s += 'Deviation: {}\n'.format(np.linalg.norm(P - np.around(P)) / 9)
        s += 'You can try raising `tolerance_cell`.'
        raise Exception(s)

    # reduce the (real) transformation matrix to the nearest integer one
    P = np.around(P)
    return P


def get_fractional_positions_from_ase_neighbor_list(atoms, neighbor_list):
    '''
    Returns the fractional positions in structure from the neighbors in the
    neighbor list.

    parameters
    ----------
    atoms : ASE Atoms object
    neighbor_list : ASE NeighborList object
    '''
    neighbor_positions = []
    fractional_positions = []

    for i in range(len(atoms)):
        lattice_site = LatticeSite_py(i, [0., 0., 0.])
        position = get_position_from_lattice_site(atoms, lattice_site)
        neighbor_positions.append(position)
        indices, offsets = neighbor_list.get_neighbors(i)
        for index, offset in zip(indices, offsets):
            lattice_site = LatticeSite_py(index, offset)
            position = get_position_from_lattice_site(atoms, lattice_site)
            neighbor_positions.append(position)
    if len(neighbor_positions) > 0:
        fractional_positions = get_scaled_positions(
            np.array(neighbor_positions),
            atoms.cell, wrap=False,
            pbc=atoms.pbc)
    return fractional_positions


def get_position_from_lattice_site(atoms, lattice_site):
    """
    Gets the corresponding position from the lattice site.

    Parameters
    ---------
    atoms : ASE atoms object
    lattice_site : icet LatticeSite object
    """
    return atoms[lattice_site.index].position + \
        np.dot(lattice_site.unitcell_offset, atoms.get_cell())


def find_lattice_site_by_position(atoms, position, tol=1e-4):
    """
    Tries to construct a lattice site equivalent from
    position in reference to the atoms object.

    atoms : ASE Atoms object
    position : x,y,z coordinate
    """

    for i, atom in enumerate(atoms):
        pos = position - atom.position
        # Direct match
        if np.linalg.norm(pos) < tol:
            return LatticeSite_py(i, np.array((0, 0, 0)))

        fractional = np.linalg.solve(atoms.cell.T, np.array(pos).T).T
        unit_cell_offset = [np.floor(round(x)) for x in fractional]
        residual = np.dot(fractional - unit_cell_offset, atoms.cell)
        if np.linalg.norm(residual) < tol:
            latNbr = LatticeSite_py(i, unit_cell_offset)
            return latNbr


def fractional_to_cartesian(atoms, frac_positions):
    """
    Turns fractional positions into cartesian positions.
    """
    return np.dot(frac_positions, atoms.cell)


def get_permutation(container, permutation):
    """
    Return the permuted version of container.
    """
    if len(permutation) != len(container):
        raise Exception
    if len(set(permutation)) != len(permutation):
        raise Exception
    return [container[s] for s in permutation]
