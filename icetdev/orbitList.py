
import numpy as np
import time
from icetdev.manybodyNeighborlist import ManybodyNeighborlist
from _icetdev import OrbitList

def __fractional_to_position(structure, row):
    """
    Return real xyz positions from fractional positions
    """
    positions = []
    for frac in row:
        position = np.dot(frac, structure.cell)
        positions.append(position)
    return positions


def __get_latNbr_permutation_matrix(structure, permutation_matrix, prune=True, verbosity=1):
    """
    Return a transformed permutation matrix with lattice neighbors instead of frac positions

    Permutation matrix is in row major format which we will keep
    """
    pm_frac = permutation_matrix.get_permutated_positions()
    pm_latNbrs = []
    for row in pm_frac:
        positions = __fractional_to_position(structure, row)
        lat_nbrs = structure.findLatticeNeighborsFromPositions(positions)
        pm_latNbrs.append(lat_nbrs)

    if prune:
        if verbosity >2:
            print("size before pruning {} ".format(len(pm_latNbrs)))
        pm_latNbrs = __prune_permutation_matrix(pm_latNbrs)
        if verbosity >2:
            print("size after pruning {} ".format(len(pm_latNbrs)))
    return pm_latNbrs


def __prune_permutation_matrix(permutation_matrix, verbosity=0):
    """
    Prunes the matrix so that the first column only contains unique elements
    """ 
    for i, rowi in enumerate(permutation_matrix):
        for j, rowj in enumerate(permutation_matrix):
            if j <= i:
                continue
            if permutation_matrix[i][0] == permutation_matrix[j][0]:
                # print("equal {} {}".format(i,j))
                permutation_matrix[j].pop()


    return permutation_matrix        
def create_orbit_list(structure, permutation_matrix, neighborlists, verbosity=3):
    """
    Create and build an orbit list from a primitive structure, a mbnl object and a permutation matrix

    structure: icet structure object (the same one used to initialize mbnl and PM
    lattice_neighbors : lattice neighbors in format given by mbnl
    permutation_matrix : icet permutation matrix object
    """

    # step1: transform permutation_matrix to be in lattice neigbhor format
    pm_lattice_neighbors =__get_latNbr_permutation_matrix(structure, permutation_matrix) 
    
    mbnl = ManybodyNeighborlist()
    t1 = time.time()
    distinct_latnbrs = mbnl.buildFromPermutationMatrix(pm_lattice_neighbors,neighborlists)
    t2 = time.time()
    if verbosity >2:
        print("time taken for building from permutation matrix {} s".format(t2-t1))
    return distinct_latnbrs

