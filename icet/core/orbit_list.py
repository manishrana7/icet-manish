from typing import List
from ase import Atoms
import numpy as np

from _icet import _OrbitList
from .local_orbit_list_generator import LocalOrbitListGenerator
from .neighbor_list import get_neighbor_lists
from .permutation_map import PermutationMap, permutation_matrix_from_atoms
from .structure import Structure
from .lattice_site import LatticeSite

from icet.io.logging import logger
logger = logger.getChild('orbit_list')


def _get_lattice_site_permutation_matrix(structure: Structure,
                                         permutation_matrix: PermutationMap,
                                         prune: bool=True):
    """
    Returns a transformed permutation matrix with lattice sites as entries
    instead of fractional coordinates.

    Parameters
    ----------
    structure
        primitive atomic structure
    permutation_matrix
        permutation matrix with fractional coordinates format entries
    prune
        if True the permutation matrix will be pruned

    Returns
    -------
    Permutation matrix in a row major order with lattice site format entries
    """
    pm_frac = permutation_matrix.get_permuted_positions()

    pm_lattice_sites = []
    for row in pm_frac:
        positions = __fractional_to_cartesian(row, structure.cell)
        lat_neighbors = []
        if np.all(structure.pbc):
            lat_neighbors = \
                structure.find_lattice_sites_by_positions(positions)
        else:
            for pos in positions:
                try:
                    lat_neighbor = \
                        structure.find_lattice_site_by_position(pos)
                except RuntimeError:
                    continue
                lat_neighbors.append(lat_neighbor)
        if len(lat_neighbors) > 0:
            pm_lattice_sites.append(lat_neighbors)
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


def __prune_permutation_matrix(permutation_matrix: List[List[LatticeSite]]):
    """
    Prunes the matrix so that the first column only contains unique elements.

    Parameters
    ----------
    permutation_matrix
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


# def create_orbit_list(atoms: Atoms, cutoffs: List[float]):
#     """
#     Builds an orbit list.

#     Parameters
#     ----------
#     atoms
#         input atomic structure
#     cutoffs
#         cutoff radii for each order

#     Returns
#     -------
#     An OrbitList object
#     """
#     max_cutoff = np.max(cutoffs)

#     # Set up a permutation matrix
#     permutation_matrix, prim_structure, _ \
#         = permutation_matrix_from_atoms(atoms, max_cutoff)

#     logger.info('Done getting permutation_matrix.')

#     # Get a list of neighbor-lists
#     neighbor_lists = get_neighbor_lists(prim_structure, cutoffs)

#     logger.info('Done getting neighbor lists.')

#     # Transform permutation_matrix to be in lattice site format
#     pm_lattice_sites \
#         = __get_lattice_site_permutation_matrix(prim_structure,
#                                                 permutation_matrix,
#                                                 prune=True)

#     msg = ['Transformation of permutation matrix to lattice neighbor']
#     msg += ['format completed.']
#     logger.info(' '.join(msg))

#     orbit_list = OrbitList(prim_structure, pm_lattice_sites, neighbor_lists)

#     logger.info('Finished construction of orbit list.')

#     return orbit_list


class OrbitList(_OrbitList):
    """
    The orbit list object has an internal list of orbits.

    An orbit has a list of equivalent sites with the restriction
    that at least one site is in the primitive cell of the
    primitive structure.

    parameters
    ----------
    atoms : ASE Atoms object
            This atoms object will be used
            to construct a primitive structure
            on which all the lattice sites in the orbits
            are based on.
    cutoffs : list of float
              cutoffs[i] is the cutoff for
              orbits with order i+2.
    """

    def __init__(self, atoms, cutoffs):
        if isinstance(atoms, Structure):
            atoms = atoms.to_atoms()
        max_cutoff = np.max(cutoffs)
        # Set up a permutation matrix
        permutation_matrix, prim_structure, _ \
            = permutation_matrix_from_atoms(atoms, max_cutoff)

        logger.info('Done getting permutation_matrix.')

        # Get a list of neighbor-lists
        neighbor_lists = get_neighbor_lists(prim_structure, cutoffs)

        logger.info('Done getting neighbor lists.')

        # Transform permutation_matrix to be in lattice site format
        pm_lattice_sites \
            = _get_lattice_site_permutation_matrix(prim_structure,
                                                   permutation_matrix,
                                                   prune=True)

        msg = ['Transformation of permutation matrix to lattice neighbor']
        msg += ['format completed.']
        logger.info(' '.join(msg))

        _OrbitList.__init__(self, prim_structure,
                            pm_lattice_sites, neighbor_lists)
        self.sort()
        logger.info('Finished construction of orbit list.')

    @property
    def primitive_structure(self):
        """
        Returns the primitive structure to
        which the lattice sites in the
        orbits are referenced to.
        """
        return self._primitive_structure.copy()

    @property
    def permutation_matrix(self):
        """
        Return icet PermutationMatrix object.
        """
        return self._permutation_matrix

    # @property
    # def orbits(self):
    #     """
    #     Return the internal list of orbits
    #     """
    #     return self._orbits

    # def __len__(self):
    #     """
    #     Lenght of an orbit list is number of orbits in orbit list.
    #     """
    #     return len(self._orbits)

    def __str__(self):
        nice_str = ''
        for i, orbit in enumerate(self.orbits):
            nice_str += "orbit {} - Multiplicity {} '\n'".format(i, len(orbit))
        return nice_str

    def get_supercell_orbit_list(self, atoms: Atoms):
        """
        Returns an orbit list for a supercell structure.

        Parameters
        ----------
        atoms
            supercell atomic structure

        Returns
        -------
        An OrbitList object

        Todo
        ----
        * Is there any reason to make this a private member
        """
        structure = Structure.from_atoms(atoms)
        log = LocalOrbitListGenerator(self, structure)

        supercell_orbit_list = log.generate_full_orbit_list()

        return supercell_orbit_list


def __fractional_to_cartesian(fractional_coordinates: List[List[float]],
                              cell: np.ndarray):
    """
    Converts cell metrics from fractional to cartesian coordinates.

    Parameters
    ----------
    fractional_coordinates
        list of fractional coordinates

    cell
        cell metric
    """
    cartesian_coordinates = [np.dot(frac, cell)
                             for frac in fractional_coordinates]
    return cartesian_coordinates
