
import numpy as np
import time
from icetdev.manybodyNeighborlist import ManybodyNeighborlist
from _icetdev import OrbitList
from icetdev.neighborlist import get_neighborlists
from icetdev.permutationMap import permutation_matrix_from_atoms

def __fractional_to_position(structure, row):
    """
    Return real xyz positions from fractional positions
    """
    positions = []
    for frac in row:
        position = np.dot(frac, structure.cell)
        positions.append(position)
    return positions


def __get_latNbr_permutation_matrix(structure, permutation_matrix, prune=True, verbosity=0):
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
        pm_latNbrs = __prune_permutation_matrix(pm_latNbrs, verbosity=0)
        if verbosity > 2:
            print("size after pruning {} ".format(len(pm_latNbrs)))
    
            
    return pm_latNbrs


def __prune_permutation_matrix(permutation_matrix, verbosity=0):
    """
    Prunes the matrix so that the first column only contains unique elements
    """ 
    rows = len(permutation_matrix)
    del_list = []
    for i  in range(len(permutation_matrix)):
        for j  in reversed(range(len(permutation_matrix))):
            if j <= i:
                continue
            if permutation_matrix[i][0] == permutation_matrix[j][0]:                
                permutation_matrix.pop(j)
                if verbosity > 2:
                    print("Removing duplicate in permutation matrix with index {1}, same as {0}".format(i,j))
    return permutation_matrix

def create_orbit_list(structure, cutoffs, verbosity=0):
    """
    Create and build an orbit list from a primitive structure, a mbnl object and a permutation matrix

    structure: icet structure object (the same one used to initialize mbnl and PM
    lattice_neighbors : lattice neighbors in format given by mbnl
    permutation_matrix : icet permutation matrix object
    """
    max_cutoff = np.max(cutoffs)
    total_time_taken = 0
    

    t0 = time.time()
    permutation_matrix, prim_structure, neighborlist = permutation_matrix_from_atoms(structure.to_atoms(), max_cutoff)    
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print("Done getting permutation_matrix. Time {0} s".format(time_taken))
    total_time_taken += time_taken

    t0 = time.time()
    neighborlists = get_neighborlists(structure=prim_structure, cutoffs=cutoffs) 
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print("Done getting neighborlists. Time {0} s".format(time_taken))




    t0 = time.time()
    #transform permutation_matrix to be in lattice neigbhor format
    pm_lattice_neighbors = __get_latNbr_permutation_matrix(prim_structure, permutation_matrix, prune=True, verbosity=verbosity)
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print("Done transforming permutation matrix into lattice neighbbor format. Time {0} s".format(time_taken))
    
    t0 = time.time()
    orbitlist = OrbitList(prim_structure, pm_lattice_neighbors, neighborlists)
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print("Done construction orbitlist. Time {0} s".format(time_taken))

    if verbosity > 3:
        print("Total time {0} s".format(total_time_taken))

    return orbitlist
   

