'''
This script checks that all atom objects in the database can have
its  clustervector computed
'''

from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms
from icetdev import permutation_map
import numpy as np
import random
from ase.db import connect


def generate_mixed_structure(atoms_prim, subelements):
    '''
    Generate a supercell structure based on the input structure and populate it
    randomly with the species specified.
    '''
    repeat = [1] * 3
    for i, pbc in enumerate(atoms_prim.pbc):
        if pbc:
            repeat[i] = 5

    atoms = atoms_prim.copy().repeat(repeat)
    for at in atoms:
        element = random.choice(subelements)
        at.symbol = element

    return atoms


def generate_clustervector_set(n, atoms_prim, subelements, clusterspace):
    '''
    Generate a set of clustervectors from clusterspace.
    '''
    clustervectors = []
    for i in range(n):
        conf = generate_mixed_structure(atoms_prim, subelements)
        conf = structure_from_atoms(conf)
        cv = clusterspace.get_clustervector(conf)
        clustervectors.append(cv)

    return clustervectors


def assert_decorrelation(matrix, tolerance=0.99):
    '''
    Confirm that the correlation between any two columns of the input matrix
    does not exceed the tolerance specified.

    Parameters
    ----------
    matrix : list of lists/NumPy array
        input matrix
    tolerance : float
        the correlation of any two columns must be lower than this value
    '''
    A = np.array(matrix)
    for i in range(len(matrix[0])):
        if i == 0:  # skip zerolets (always one)
            continue
        col_i = matrix[:, i]
        for j in range(len(matrix[0])):
            if j <= i:
                continue
            col_j = A[:, j]
            corr = np.dot(col_i, col_j)
            corr /= np.linalg.norm(col_i) * np.linalg.norm(col_j)
            msg = 'columns {} and {} correlated <i|j>= {}'.format(i, j, corr)
            assert corr < tolerance, msg


print('')
db = connect('structures_for_testing.db')
subelements = ['H', 'He', 'Pb']
for row in db.select():
    atoms_row = row.toatoms()
    atoms_tag = row.tag
    cutoffs = [1.4] * 3
    if len(atoms_row) == 0:
        continue
    if atoms_row.get_pbc().all():
        atoms_row.wrap()
        print(' structure: {}'.format(row.tag))
        clusterspace = create_clusterspace(atoms_row, cutoffs, subelements)
        if not atoms_row.get_pbc().all():
            permutation_map.__vacuum_on_non_pbc(atoms_row)
        cvs = generate_clustervector_set(5, atoms_row,
                                         subelements, clusterspace)
