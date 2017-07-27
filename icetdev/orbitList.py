
import numpy as np


def __fractional_to_position(structure, row):
    """
    Return real xyz positions from fractional positions
    """
    positions = []
    for frac in row:
        position = np.dot(frac, structure.cell)
        positions.append(position)
    return positions


def __get_latNbr_permutation_matrix(structure, permutation_matrix):
    """
    Return a transformed permutation matrix with latnbrs instead of frac positions

    Permutation matrix is in row major format which we will keep
    """
    pm_frac = permutation_matrix.get_permutated_positions()
    pm_latNbrs = []
    for row in pm_frac:
        positions = __fractional_to_position(structure, row)
        lat_nbrs = structure.findLatticeNeighborsFromPositions(positions)
        pm_latNbrs.append(lat_nbrs)

    return pm_latNbrs


def create_orbit_list(structure, mbnl, permutation_matrix):
    """
    Create and build an orbit list from a primitive structure, a mbnl object and a permutation matrix

    structure: icet structure object (the same one used to initialize mbnl and PM
    mbnl : icet manybodyNeighborlist object
    permutation_matrix : icet permutation matrix object
    """

    # step1: transform permutation_matrix to be in lattice neigbhor format
    pm_latticeNeighbors =__get_latNbr_permutation_matrix(structure, mbnl) 
    