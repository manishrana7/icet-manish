import time
import numpy as np
from _icet import OrbitList
from .local_orbit_list_generator import LocalOrbitListGenerator
from .neighbor_list import get_neighbor_lists
from .permutation_map import permutation_matrix_from_atoms
from .structure import Structure

from icet.io.logging import logger
logger = logger.getChild('orbit_list')


def __fractional_to_cartesian(fractional_coordinates, cell):
    """
    Convert from fractional to Cartesian coordinates.

    Parameters
    ----------
    fractional_coordinates : list of 3d-vectors
        list fractional coordinates

    cell : 3x3 matrix
        cell metric
    """
    cartesian_coordinates = [np.dot(frac, cell)
                             for frac in fractional_coordinates]
    return cartesian_coordinates


def __get_lattice_site_permutation_matrix(structure, permutation_matrix,
                                          prune=True):
    """
    Return a transformed permutation matrix with lattice sites instead of
    fractional coordinates.

    Parameters
    ----------
    structure : icet Structure object
        primitive atomic structure

    permutation_matrix : icet PermutatioMap object
        permutation matrix

    prune : bool
        if True prune the permutation matrix. Default to True

    Permutation matrix is in row major format which we will keep
    """
    pm_frac = permutation_matrix.get_permuted_positions()

    pm_lattice_sites = []
    for row in pm_frac:
        positions = __fractional_to_cartesian(row, structure.cell)
        lat_nbrs = []
        if np.all(structure.pbc):
            lat_nbrs = structure.find_lattice_sites_by_positions(positions)
        else:
            for pos in positions:
                try:
                    lat_nbr = structure.find_lattice_site_by_position(pos)
                    lat_nbrs.append(lat_nbr)
                except:  # NOQA
                    continue
        if len(lat_nbrs) > 0:
            pm_lattice_sites.append(lat_nbrs)
        else:
            logger.warning('Unable to transform any element in a column of the'
                           ' fractional permutation matrix to lattice site')
    if prune:
        logger.debug('Size of columns of the permutation matrix before'
                     ' pruning {}'.format(len(pm_lattice_sites)))

        pm_lattice_sites = __prune_permutation_matrix(pm_lattice_sites)

        logger.debug('Size of columns of the permutation matrix after'
                     ' pruning {}'.format(len(pm_lattice_sites)))

    return pm_lattice_sites


def __prune_permutation_matrix(permutation_matrix):
    """
    Prunes the matrix so that the first column only contains unique elements.

    Parameters
    ----------
    permutation_matrix : matrix
        permutation matrix with LatticeSite type entries
    """
    for i in range(len(permutation_matrix)):
        for j in reversed(range(len(permutation_matrix))):
            if j <= i:
                continue
            if permutation_matrix[i][0] == permutation_matrix[j][0]:
                permutation_matrix.pop(j)
                msg = ['Removing duplicate in permutation matrix']
                msg += ['i: {} j: {}'.format(i, j)]
                logger.debug(' '.join(msg))

    return permutation_matrix


def _get_supercell_orbit_list(self, atoms):
    """
    Returns an orbit list for a supercell structure

    Parameters
    ----------
    atoms: ASE Atoms object
        supercell atomic structure
    """
    structure = Structure.from_atoms(atoms)
    log = LocalOrbitListGenerator(self, structure)

    supercell_orbit_list = log.generate_full_orbit_list()

    return supercell_orbit_list


OrbitList.get_supercell_orbit_list = _get_supercell_orbit_list


def create_orbit_list(structure, cutoffs):
    '''
    Build an orbit list.

    Parameters
    ----------
    structure: icet Structure or ASE Atoms object (bioptional)
        input configuration used to initialize mbnl and permutation matrix
    cutoffs : list of float
        cutoff radii for each order

    Returns
    -------
    OrbitList object
    '''
    if isinstance(structure, Structure):
        atoms = Structure.to_atoms(structure)
    else:
        atoms = structure.copy()

    max_cutoff = np.max(cutoffs)
    total_time_spent = 0

    t0 = time.time()
    permutation_matrix, prim_structure, neighbor_list \
        = permutation_matrix_from_atoms(atoms, max_cutoff)
    time_spent = time.time() - t0
    total_time_spent += time_spent

    msg = 'Done getting permutation_matrix (time: {:.6f}s)'.format(time_spent)
    logger.info(msg)

    t0 = time.time()
    neighbor_lists = get_neighbor_lists(prim_structure, cutoffs=cutoffs)
    t1 = time.time()
    time_spent = t1 - t0
    total_time_spent += time_spent

    logger.info('Done getting neighbor lists.'
                ' (time: {:.6f}s)'.format(time_spent))

    t0 = time.time()
    # transform permutation_matrix to be in lattice site format
    pm_lattice_sites \
        = __get_lattice_site_permutation_matrix(prim_structure,
                                                permutation_matrix,
                                                prune=True)
    t1 = time.time()
    time_spent = t1 - t0
    total_time_spent += time_spent

    msg = ['Transformation of permutation matrix to lattice neighbor']
    msg += ['format completed (time: {:.6f}s)'.format(time_spent)]
    logger.info(' '.join(msg))

    t0 = time.time()
    orbit_list = OrbitList(prim_structure, pm_lattice_sites, neighbor_lists, bothways)
    t1 = time.time()
    time_spent = t1 - t0
    total_time_spent += time_spent

    logger.info('Finished construction of orbit list.'
                ' (time: {:.6f}s)'.format(time_spent))

    logger.info('Total time: {:.6f}s'.format(total_time_spent))

    return orbit_list
