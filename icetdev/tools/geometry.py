import numpy as np
from icetdev.lattice_site import LatticeSite
from icetdev import Structure
import math

from collections import namedtuple


def get_scaled_positions(positions, cell, wrap=True, pbc=[True, True, True]):
    """Get positions relative to unit cell.

    If wrap is True, positions outside the unit cell will be wrapped into
    the cell in those directions with periodic boundary conditions
    so that the scaled coordinates are between zero and one.
    """

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
    """
    Get lattice neighbor from position

    This is the Python version of structure.findLatticeSiteFromPosition(position)

    It is slower but kept as help for debugging and if further development is needed
    """

    fractional = np.linalg.solve(structure.cell.T, np.array(position).T).T
    unit_cell_offset = [int(round(x)) for x in fractional]

    remainder = np.dot(fractional - unit_cell_offset, structure.cell)
    try:
        index = structure.find_index_of_position(remainder, 1e-3)
    except:
        print("error did not find index with pos: {}".format(remainder))
        print("position in structure are:")
        print(structure.positions)
        exit(1)

    latNbr = LatticeSite(index, unit_cell_offset)
    return latNbr


# def transform_cell_to_cell(atoms, atoms_template):
#     '''
#     Transform atoms_transform to look like a simple repeat of
#     atoms_template.
#     '''

#     atoms = atoms.copy()

#     cell_transform = atoms.cell
#     cell_template = atoms_template.cell

#     #get rotation matrix to rotate atoms into atoms_template
#     rotation_matrix = np.linalg.solve(cell_transform, cell_template)/np.linalg.norm(cell_transform )
#     # rotation_matrix = np.linalg.solve(cell_transform,cell_transform)/np.linalg.norm(cell_transform )

#     atoms.cell = np.dot(cell_transform, cell_template)/np.linalg.norm(cell_template)

#     for atom in atoms:
#         atom.position = np.dot(rotation_matrix, atom.position)

#     return atoms

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
        offset = tuple(np.floor(np.round(pos, decimals=5)).astype(int))
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


def get_permutation_matrix(input_configuration,
                           reference_structure,
                           tolerance_cell=0.05,
                           ):
    '''
    Computes and returns the permutation matrix that takes the reference cell to the input cell,
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
    smart_offsets: List of named tuple
        each named tuple is one prmi

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
            smart_offsets.append([offset] * size_of_primitive)
        return smart_offsets

    # Sanity checks
    if len(unique_offsets) * size_of_primitive < size_of_supercell:
        raise Exception("Undefined behaviour in function get_smart_offsets")


    # Initialize map
    supercell_maps = {}
    for i, p in enumerate(atoms.positions):
        supercell_maps[i] = []

    """
    supercell[i] = offset means that supercell atom `i` was mapped by some primitive atom
    when the primitive is transleted by offset

    """

    offset_mapping_list = []
    primitive_map = {}
    for i in range(len(atoms_prim)):
        primitive_map[i] = []
    for offset in unique_offsets:
        positions = get_offset_positions(atoms_prim, offset)
        for prim_index, pos in enumerate(positions):

            # Difference between the offset of the primitive position and supercell positions
            pos_diff = np.floor(np.round(atoms.positions - pos, decimals=7))
            matched_indices = get_indices_with_zero_component(pos_diff)
            for i in matched_indices:
                mapped_entry = [prim_index, i, offset]
                offset_mapping_list.append(mapped_entry)
                supercell_maps[i].append([offset, prim_index])
                primitive_map[prim_index].append([offset, i])

    OffsetMap = namedtuple("OffsetMap", "prim_indices, super_indices, offsets")

    OneOffsetMat = namedtuple(
        "OneOffsetMat", "prim_index, super_index, offset")


    """
    Construct smart offsets so that
    1: the list of offset are have a size of len(supercell)/len(primitive)
    2: The different offsets will translate the primitive so that
       each interatomic distance will be kept (under the supercells cell)
    3: each offset in the list of offset will map to an atom in the supercell which
       correspond to a lattice site that has the same index is the translated atom

    Construct all interatomic distances of the primitive cell.

    When constructing the elements in list of offsets, check that the interatomic distances are kept

    """
    number_of_offsets = len(atoms) // len(atoms_prim)
    list_of_offset = []
    list_of_tuple_offset = []
    print("number of offsets", number_of_offsets)
    for kk in range(number_of_offsets):
        current_offsets = []
        currentOffsetsTuple = OffsetMap([], [], [])

        for offset_index in reversed(offset_mapping_list):
            next_offset = OneOffsetMat(
                offset_index[0], offset_index[1], offset_index[2])

            if is_compatible_new_offset(currentOffsetsTuple, next_offset, atoms, atoms_prim):
                currentOffsetsTuple.prim_indices.append(offset_index[0])
                currentOffsetsTuple.super_indices.append(offset_index[1])
                currentOffsetsTuple.offsets.append(offset_index[2])
                offset_mapping_list.pop(
                    offset_mapping_list.index(offset_index))

            if len(currentOffsetsTuple.super_indices) == len(atoms_prim):
                break
        assert len(currentOffsetsTuple.prim_indices) == len(atoms_prim)
        current_offsets = []
        for i, j, offset in zip(currentOffsetsTuple.prim_indices, currentOffsetsTuple.super_indices, currentOffsetsTuple.offsets):
            current_offsets.append([i, j, offset])
        list_of_tuple_offset.append(currentOffsetsTuple)
        list_of_offset.append(current_offsets)
    return(list_of_tuple_offset)


def is_compatible_new_offset(current_offsets, next_offsets, atoms, atoms_prim, tol=1e-5):
    """
    Tests if next offsets are compatible to be added to current offsets
    """
    if len(current_offsets) == 0:
        return True
    # See that this primitive index has not been mapped in current
    if next_offsets.prim_index in current_offsets.prim_indices:
        return False
    # See that this supercell index has not been mapped in current
    if next_offsets.super_index in current_offsets.super_indices:
        return False

    for i, j in zip(current_offsets.prim_indices, current_offsets.super_indices):
        dist_prim = atoms_prim.get_distance(
            i, next_offsets.prim_index, mic=True)
        dist_super = atoms.get_distance(j, next_offsets.super_index, mic=True)
        if np.abs(dist_prim - dist_super) > tol:
            return False
    return True


def get_unitcell_offsets_from_positions(positions, cell):
    """
    Return the unitcell that the positions resides in.

    Parameters
    ----------
    positions : list/array of xyz positions

    cell : numpy 3x3 array
        the unitcell to relate the positions to
    """
    fractional_positions = get_scaled_positions(positions, cell, wrap=False)
    unit_cells = []
    for frac_pos in fractional_positions:
        offset = tuple(np.floor(np.round(frac_pos, decimals=5)))
        unit_cells.append(offset)
    return unit_cells


def get_indices_with_zero_component(array, tolerance=1e-4):
    """
    Returns the indices of the array where the component is zero within a tolerance
    """
    indices = []
    for i, a in enumerate(array):
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

    return np.array(atoms.positions + np.dot(offset, atoms.cell))
