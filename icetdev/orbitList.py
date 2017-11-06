
import numpy as np
import time
from _icetdev import OrbitList
from icetdev.neighborlist import get_neighborlists
from icetdev.permutationMap import permutation_matrix_from_atoms


def __fractional_to_cartesian(fractional_coordinates, cell):
    '''
    Convert from fractional to Cartesian coordinates.

    Parameters
    ----------
    fractional_coordinates : list of 3d-vectors
        list fractional coordinates
    cell : 3x3 matrix
        cell metric
    '''
    cartesian_coordinates = [np.dot(frac, cell)
                             for frac in fractional_coordinates]
    return cartesian_coordinates


def __get_latNbr_permutation_matrix(structure, permutation_matrix,
                                    prune=True, verbosity=0):
    '''
    Return a transformed permutation matrix with lattice neighbors instead of
    fractional coordinates.

    Permutation matrix is in row major format which we will keep
    '''
    pm_frac = permutation_matrix.get_permutated_positions()

    pm_latNbrs = []
    for row in pm_frac:
        positions = __fractional_to_cartesian(row, structure.cell)
        lat_nbrs = []
        if np.all(structure.pbc):
            lat_nbrs = structure.find_lattice_neighbors_by_positions(positions)
        else:
            for pos in positions:
                try:
                    lat_nbr = structure.find_lattice_neighbor_by_position(pos)
                    lat_nbrs.append(lat_nbr)
                except:
                    continue
        if len(lat_nbrs) > 0:
            pm_latNbrs.append(lat_nbrs)
        else:
            print('lat nbrs are zero')
    if prune:
        if verbosity > 2:
            print('size before pruning {} '.format(len(pm_latNbrs)))
        pm_latNbrs = __prune_permutation_matrix(pm_latNbrs, verbosity=0)
        if verbosity > 2:
            print('size after pruning {} '.format(len(pm_latNbrs)))

    return pm_latNbrs


def __prune_permutation_matrix(permutation_matrix, verbosity=0):
    '''
    Prunes the matrix so that the first column only contains unique elements
    '''
    for i in range(len(permutation_matrix)):
        for j in reversed(range(len(permutation_matrix))):
            if j <= i:
                continue
            if permutation_matrix[i][0] == permutation_matrix[j][0]:
                permutation_matrix.pop(j)
                if verbosity > 2:
                    msg = ['Removing duplicate in permutation matrix']
                    msg += ['i: {} j: {}'.format(i, j)]
                    print(' '.join(msg))
    return permutation_matrix


def create_orbit_list(structure, cutoffs, verbosity=0):
    '''
    Build an orbit list from a primitive structure, a many-body neighbor list
    (mbnl) object, and a permutation matrix.

    Parameters
    ----------
    structure: icet structure object
        input configuration used to initialize mbnl and permutation matrix
    lattice_neighbors : mbnl object
        lattice neighbors
    permutation_matrix : icet permutation matrix object
        permutation matrix
    '''
    max_cutoff = np.max(cutoffs)
    total_time_taken = 0

    t0 = time.time()
    permutation_matrix, prim_structure, neighborlist \
        = permutation_matrix_from_atoms(structure.to_atoms(), max_cutoff)
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print('Done getting permutation_matrix. Time {} s'.format(time_taken))
    total_time_taken += time_taken

    t0 = time.time()
    neighborlists = get_neighborlists(structure=prim_structure,
                                      cutoffs=cutoffs)
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print('Done getting neighborlists. Time {} s'.format(time_taken))

    t0 = time.time()
    # transform permutation_matrix to be in lattice neigbhor format
    pm_lattice_neighbors \
        = __get_latNbr_permutation_matrix(prim_structure, permutation_matrix,
                                          prune=True, verbosity=verbosity)
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        msg = ['Transformation of permutation matrix to lattice neighbor']
        msg += ['format completed (time: {} s)'.format(time_taken)]
        print(' '.join(msg))

    t0 = time.time()
    orbitlist = OrbitList(prim_structure, pm_lattice_neighbors, neighborlists)
    t1 = time.time()
    time_taken = t1 - t0
    total_time_taken += time_taken

    if verbosity > 3:
        print('Done construction orbitlist. Time {0} s'.format(time_taken))

    if verbosity > 3:
        print('Total time {0} s'.format(total_time_taken))

    return orbitlist
