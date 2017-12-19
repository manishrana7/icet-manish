import numpy as np
from ase import Atoms
import spglib

from icetdev.core.lattice_site import LatticeSite


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

<<<<<<< 1c0f42508bd7362df7e598120a7d00c215181349
=======
    It is slower but kept as help for debugging and if further development is needed
    """

>>>>>>> DEV #55: staking out the bug fix in python
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

<<<<<<< 1c0f42508bd7362df7e598120a7d00c215181349
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
=======
# def transform_cell_to_cell(atoms, atoms_template):
#     '''
#     Transform atoms_transform to look like a simple repeat of
#     atoms_template.
#     '''
>>>>>>> DEV #55: staking out the bug fix in python

    if len(vacuum_axis) > 0:
        atoms.center(30, axis=vacuum_axis)
    atoms.wrap()

    return atoms


def get_primitive_structure(atoms):
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
        atoms_with_vacuum, to_primitive=True, no_idealize=True)
    scaled_positions = [np.round(pos, 12) for pos in scaled_positions]
    atoms_prim = Atoms(scaled_positions=scaled_positions,
                       numbers=numbers, cell=lattice, pbc=atoms.pbc)
    atoms_prim.wrap()
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

def required_offsets_to_map_supercell(supercell, atoms_prim):
    '''
    Calculates the minimum number of offsets 
    of atoms prim needed to completely cover the atoms object

    Parameters
    ----------
    supercell: ASE atoms object
        The supercell
    atoms_prim: ASE atoms object
        The primitive cell

    Returns
    ------
    required_offsets: List of lists
        A minimum set of offsets of the primitive 
        needed to cover the supercell
    '''
    # Get fractional coordinates of supercell positions given in primitive cell
    fractional_positions = get_scaled_positions(
        supercell.positions, atoms_prim.cell, wrap=False)

    offsets = []

    for pos in fractional_positions:
        offset = tuple(np.floor(pos).astype(int))
        offsets.append(offset)

    required_offsets = list(set(offsets))

    return required_offsets


def transform_cell_to_cell(atoms, atoms_template, tolerance=1e-3):
    '''
    Transform atoms_transform to look like a simple repeat of
    atoms_template.
    '''

    atoms = atoms.copy()
    atoms_template = atoms_template.copy()
    cell_transform = atoms.cell
    cell_template = atoms_template.cell
    atoms.wrap()
    atoms_template.wrap()

    # get fractional coordinates of supercell positions relative primitive cell
    fractional_positions = get_scaled_positions(
        atoms.positions, atoms_template.cell, wrap=False)

    # print(atoms.positions)

    offsets = []

    for pos in fractional_positions:
        offset = tuple(np.floor(pos).astype(int))
        offsets.append(offset)

    unique_offsets = list(set(offsets))
    max_offset = list(max(unique_offsets))
    for i in range(3):
        if max_offset[i] < 0:
            max_offset[i] = 0
        max_offset[i] += 1
    print(max_offset)

    atoms_new = atoms_template.copy().repeat(max_offset)

    scaled_positions = atoms_new.get_scaled_positions()

    print(atoms_new.cell)
    supercell_fractional = get_scaled_positions(
        atoms.positions,  atoms_new.cell, wrap=True)

    for i, atom in enumerate(atoms_new):
        for j in range(len(supercell_fractional)):
            if np.linalg.norm(scaled_positions[i] - supercell_fractional[j]) < tolerance:
                atom.symbol = atoms[j].symbol

    return atoms_new

<<<<<<< 1c0f42508bd7362df7e598120a7d00c215181349
=======

>>>>>>> DEV #55: staking out the bug fix in python
def get_permutation_matrix(input_configuration,
                           reference_structure,
                           tolerance_cell=0.05,
                           ):
    '''
    Computes and returns the permutation matrix that takes the reference cell to the input cell,
    i.e. permutation_matrix * reference_cell = input_cell
    '''

    input_cell = input_configuration.cell
    reference_cell = reference_structure.cell

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


def get_smart_offsets(atoms, atoms_prim):
    '''
     Returns the maps (note plural) that maps each basis atom to the supercell once and only once
     Each basis will get its own offset. 

    Parameters
    ----------
    atoms: ASE atoms object
        The supercell
    atoms_prim: ASE atoms object
        The primitive cell

    Returns
    ------
    smart_offsets: List of list
        each inner list contains a subset of the mapping
        needed to map the entire supercell from the primitive cell.

    Raises
        Exception 
            if the algorithm doesn't find any mapping to the supercelll
    '''

    smart_offsets = []

    size_of_primitive = len(atoms_prim)
    size_of_supercell = len(atoms)
    expected_number_of_mappings = size_of_supercell // size_of_primitive
    assert np.abs(expected_number_of_mappings - size_of_supercell /
                  size_of_primitive) < 0.01, "Can not translate cell with {} atoms to a supercell with {} atoms".format(size_of_primitive, size_of_supercell)

    unique_offsets = required_offsets_to_map_supercell(atoms, atoms_prim)

    # Check if easy solution is possibly
    if len(unique_offsets) == expected_number_of_mappings:
        for offset in unique_offsets:
            smart_offsets.append( [offset] * size_of_primitive )
        return smart_offsets               
    
    # Sanity checks
    if len(unique_offsets) * size_of_primitive < size_of_supercell:
        raise Exception("Undefined behaviour in function get_smart_offsets")


    mapped_supercell_atoms = []
    
    # Initialize map
    supercell_maps = {}
    for i,p in enumerate(atoms.positions):
        supercell_maps[i] = []


    for offset in unique_offsets:
        positions = get_offset_positions(atoms_prim, offset)
        for pos in positions:
            matched_indices = get_indices_with_zero_component(atoms.positions - pos)
            for i in matched_indices:
                supercell_maps[i].append(offset)
                

    for key in supercell_maps.keys():
        print(key, end= ': ')
        for offsets in supercell_maps[key]:
            print(offsets, end=', ')
        print()            
    return smart_offsets 



def get_indices_with_zero_component(array, tolerance=1e-4):
    """ 
    Returns the indices of the array where the component is zero within a tolerance
    """
    indices = []
    for i,a in enumerate(array):
        if np.linalg.norm(a) < tolerance:
            indices.append(i)
    return indices            



def get_offset_positions(atoms, offset):
    ''' 
    Get the offset positions
    
    parameters
    ---------
    atoms : ASE Atoms object
        atoms object from which the positions are to be taken
        and the offsets should be related to its unitcell
    offset: list of int 
        Offsets of the unitcell vectors in [O_x, O_y, O_z]
     '''
    return atoms.positions + np.dot(offset, atoms.cell)