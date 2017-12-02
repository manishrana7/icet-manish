import numpy as np
from icetdev.lattice_site import LatticeSite
import math

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
        index = structure.find_index_of_position(remainder)
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

    Paramaters
    ----------
    supercell: ASE atoms object
        The supercell
    atoms_prim: ASE atombs object
        The primitive cell

    Returns
    ------
    required_offsets: List of lists
        A minimum set of offsets of the primitive 
        needed to cover the supercell
    '''
    # Get fractional coordinates of supercell positions given in primitive cell
    fractional_positions = get_scaled_positions(supercell.positions, atoms_prim.cell, wrap=False)

    offsets = []

    for pos in fractional_positions:
        offset = tuple(np.floor(pos).astype(int))
        offsets.append(offset)

    required_offsets = list(set(offsets))
    
    return required_offsets

def transform_cell_to_cell(atoms, atoms_template, tolerance = 1e-3):
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
    fractional_positions = get_scaled_positions(atoms.positions, atoms_template.cell, wrap=False)

    # print(atoms.positions)
    
    offsets = []

    for pos in fractional_positions:
        offset = tuple(np.floor(pos).astype(int))
        offsets.append(offset)


    unique_offsets = list(set(offsets))
    max_offset = list(max(unique_offsets))
    for i in range(3):
        if max_offset[i] < 0 :
            max_offset[i] = 0
        max_offset[i] +=1
    print(max_offset)    
    
    atoms_new = atoms_template.copy().repeat(max_offset)
    
    scaled_positions = atoms_new.get_scaled_positions()

    print(atoms_new.cell)
    supercell_fractional = get_scaled_positions(atoms.positions,  atoms_new.cell, wrap=True)

    for i, atom in enumerate(atoms_new):
        for j in range( len(supercell_fractional)):
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