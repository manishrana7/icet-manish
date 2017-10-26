"""
This script checks the column correlation for a set of clustervectors,
it will be asserted that no columns are not highly correlated
"""

from icetdev import clusterspace
from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms
from ase.build import bulk, make_supercell
from icetdev import permutationMap
import numpy as np
import random
from ase.db import connect


def generateRandomStructure(atoms_prim, subelements):
    """
    Generate a random structure with atoms_prim as a base
    and fill it randomly with elements in subelements
    """

    repeat = [1]*3
    for i, pbc in enumerate(atoms_prim.pbc):
        if pbc:
            repeat[i] = 5

    
    atoms = atoms_prim.copy().repeat(repeat)
    
    for at in atoms:
        element = random.choice(subelements)
        at.symbol = element

    if len(atoms) < len(atoms_prim):
        raise Exception("Error: generate random structure produced a smaller system")
    return atoms


def generateCVSet(n, atoms_prim, subelements, clusterspace):
    """
    generate a set of clustervectors from clusterspace
    """
    clustervectors = []
    for i in range(n):
        conf = generateRandomStructure(atoms_prim, subelements)        
            
        conf = structure_from_atoms(conf)
        cv = clusterspace.get_clustervector(conf)
        clustervectors.append(cv)

    return clustervectors


def getColumnCorrelation(i, j, cv_matrix):
    """
    Returns the correlation between column i and j

    cv_matrix: numpy matrix
    """
    col_i = cv_matrix[:, i]
    col_j = cv_matrix[:, j]
    
    
    corr = np.dot(col_i, col_j) / \
        (np.linalg.norm(col_i) * np.linalg.norm(col_j))
          
    return corr


def assertNoCorrelation(cvs, tol=0.99):
    """
    check that no column in cvs are above tolerance
    """
    cvs_matrix = np.array(cvs)
    for i in range(len(cvs[0])):
        if i == 0: #dont loop over zerolet since this is always ones
            continue
        for j in range(len(cvs[0])):
            if j <= i:
                continue
            corr = getColumnCorrelation(i, j, cvs_matrix)
            assert corr < tol, "columns {} and {} were correletated with {}".format(
                i, j, corr)


db = connect("structures_for_testing.db")

subelements = ["H", "He", "Pb"]
for row in db.select():
    atoms_row = row.toatoms()
    atoms_tag = row.tag
    cutoffs = [1.4] * 3
    if len(atoms_row) == 0:
        continue
    if atoms_row.get_pbc().all() == True:
        atoms_row.wrap()
    if True: # atoms_row.get_pbc().all() == True:
        print("Testing finding cv for structure: {} with cutoffs {}".format(atoms_tag, cutoffs))
        
        clusterspace = create_clusterspace(
            subelements, cutoffs, atoms=atoms_row)
        if not atoms_row.get_pbc().all() == True:
            permutationMap.__vacuum_on_non_pbc(atoms_row)

        cvs = generateCVSet(5, atoms_row, subelements, clusterspace)

