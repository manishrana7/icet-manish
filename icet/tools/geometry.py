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
            return LatticeSite_py(i, np.array((0., 0., 0.)))

        fractional = np.linalg.solve(atoms.cell.T, np.array(pos).T).T
        unit_cell_offset = [np.floor(round(x)) for x in fractional]
        residual = np.dot(fractional - unit_cell_offset, atoms.cell)
        if np.linalg.norm(residual) < tol:
            latNbr = LatticeSite_py(i, unit_cell_offset)
            return latNbr

    # found nothing, raise error
    raise RuntimeError("Did not find site in find_lattice_site_by_position")


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
        raise RuntimeError("Containter and permutation"
                           " not of same size {} != {}".format(
                               len(container), len(permutation)))
    if len(set(permutation)) != len(permutation):
        raise Exception
    return [container[s] for s in permutation]


def find_permutation(target, permutated):
    """
    Returns the permutation vector that takes
    permutated to target

    parameters
    ----------
    target : some container
        container should allow .index and the
    containers elements should contain objects
    with __eq__ method
    permutated : some container
        container should allow .index and the
    containers elements should contain objects
    with __eq__ method
    """
    permutation = []
    for element in target:
        index = permutated.index(element)
        permutation.append(index)
    return permutation
