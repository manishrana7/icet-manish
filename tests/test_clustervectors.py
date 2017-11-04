'''
This script checks that all atom objects in the database can have
its  clustervector computed
'''

from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms
from icetdev import permutationMap
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


def get_column_correlation(i, j, matrix):
    '''
    Returns the correlation between two columns i and j of the input matrix.

    Parameters
    ----------
    i : int
        index of first column
    j : int
        index of second column
    matrix: NumPy array
        input matrix
    '''
    col_i = matrix[:, i]
    col_j = matrix[:, j]
    corr = np.dot(col_i, col_j)
    corr /= np.linalg.norm(col_i) * np.linalg.norm(col_j)
    return corr


def assert_decorrelation(matrix, tol=0.99):
    '''
    Confirm that the correlation between any two columns of the input matrix
    does not exceed the tolerance specified.
    '''
    A = np.array(matrix)
    for i in range(len(matrix[0])):
        if i == 0:  # skip zerolets (always one)
            continue
        for j in range(len(matrix[0])):
            if j <= i:
                continue
            corr = get_column_correlation(i, j, A)
            msg = 'columns {} and {} correlated <i|j>= {}'.format(i, j, corr)
            assert corr < tol, msg


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
            permutationMap.__vacuum_on_non_pbc(atoms_row)
        cvs = generate_clustervector_set(5, atoms_row,
                                         subelements, clusterspace)