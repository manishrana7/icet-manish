"""
This module provides the OrbitList class.
"""

from typing import List, OrderedDict, Dict

import numpy as np

from _icet import _OrbitList
from ase import Atoms
from icet.core.local_orbit_list_generator import LocalOrbitListGenerator
from icet.core.neighbor_list import get_neighbor_lists
from icet.core.matrix_of_equivalent_positions import \
    _get_lattice_site_matrix_of_equivalent_positions, \
    matrix_of_equivalent_positions_from_structure
from icet.core.structure import Structure
from icet.input_output.logging_tools import logger
from icet.tools.geometry import atomic_number_to_chemical_symbol

logger = logger.getChild('orbit_list')


class OrbitList(_OrbitList):
    """
    The orbit list object handles an internal list of orbits.

    An orbit has a list of equivalent sites with the restriction
    that at least one site is in the cell of the primitive structure.

    Parameters
    ----------
    structure : Atoms
        This structure will be used to construct a primitive
        structure on which all the lattice sites in the orbits
        are based.
    cutoffs : list of float
        the i-th element of this list is the cutoff for orbits with
        order i+2
    symprec : float
        tolerance imposed when analyzing the symmetry using spglib
    position_tolerance : float
        tolerance applied when comparing positions in Cartesian coordinates
    fractional_position_tolerance : float
        tolerance applied when comparing positions in fractional coordinates
    """

    def __init__(self,
                 structure: Atoms,
                 cutoffs: List[float],
                 symprec: float,
                 position_tolerance: float,
                 fractional_position_tolerance: float) -> None:
        max_cutoff = np.max(cutoffs)
        # Set up a permutation matrix
        matrix_of_equivalent_positions, prim_structure, _ \
            = matrix_of_equivalent_positions_from_structure(structure=structure,
                                                            cutoff=max_cutoff,
                                                            position_tolerance=position_tolerance,
                                                            find_primitive=False,
                                                            symprec=symprec)

        logger.info('Done getting matrix_of_equivalent_positions.')

        # Get a list of neighbor-lists
        neighbor_lists = get_neighbor_lists(structure=prim_structure, cutoffs=cutoffs,
                                            position_tolerance=position_tolerance)

        logger.info('Done getting neighbor lists.')

        # Transform matrix_of_equivalent_positions to be in lattice site format
        pm_lattice_sites = _get_lattice_site_matrix_of_equivalent_positions(
            structure=prim_structure,
            matrix_of_equivalent_positions=matrix_of_equivalent_positions,
            fractional_position_tolerance=fractional_position_tolerance,
            prune=True)

        logger.info('Transformation of matrix of equivalent positions'
                    ' to lattice neighbor format completed.')

        _OrbitList.__init__(self,
                            structure=prim_structure,
                            matrix_of_equivalent_sites=pm_lattice_sites,
                            neighbor_lists=neighbor_lists,
                            position_tolerance=position_tolerance)
        logger.info('Finished construction of orbit list.')

    @property
    def primitive_structure(self):
        """
        Returns the primitive structure to which the lattice sites in
        the orbits are referenced to.
        """
        return self._primitive_structure.copy()

    def __str__(self):
        """String representation."""
        s = []
        s += ['Number of orbits: {}'.format(len(self))]
        for k, orbit in enumerate(self.orbits):
            c = self.orbits[k].representative_cluster.__str__()
            s += [f'orbit: {k:3}   order: {orbit.order:3}'
                  f'   multiplicity: {len(orbit):3}   representative_cluster: {c}']
        return '\n'.join(s)

    def get_supercell_orbit_list(self,
                                 structure: Atoms,
                                 fractional_position_tolerance: float):
        """
        Returns an orbit list for a supercell structure.

        Parameters
        ----------
        structure
            supercell atomic structure
        fractional_position_tolerance : float
            tolerance applied when comparing positions in fractional coordinates
        """
        lolg = LocalOrbitListGenerator(
            self,
            structure=Structure.from_atoms(structure),
            fractional_position_tolerance=fractional_position_tolerance)
        supercell_orbit_list = lolg.generate_full_orbit_list()
        return supercell_orbit_list

    def remove_inactive_orbits(self,
                               allowed_species: List[List[str]]) -> None:
        """ Removes orbits with inactive sites.

        Parameters
        ----------
        allowed_species
            the list of allowed species on each site in the primitive
            structure
        """
        prim_structure = self.get_primitive_structure()
        number_of_allowed_species = [len(sym) for sym in allowed_species]
        prim_structure.set_number_of_allowed_species(number_of_allowed_species)
        self._remove_inactive_orbits(prim_structure)

    def get_cluster_counts(self,
                           structure: Atoms,
                           fractional_position_tolerance,
                           orbit_indices: List[int] = None) -> OrderedDict[int, Dict]:
        """
        Counts all clusters in a structure by finding their local orbit list.

        Parameters
        ----------
        structure
            Supercell (commensurate with the structure this orbit list is
            based on) to count clusters for.
        fractional_position_tolerance : float
            tolerance applied when comparing positions in fractional coordinates
        orbit_indices
            Indices of orbits for which counts are desired. If None, all
            orbits will be counted.
        
        Returns
        -------
        Ordered dict, the keys of which are orbit indices and the values
        cluster counts. The latter are themselves dicts, with tuples
        of chemical symbols as keys and the number of such clusters
        as values.
        """
        supercell_orbit_list = self.get_supercell_orbit_list(
            structure=structure,
            fractional_position_tolerance=fractional_position_tolerance)

        # Collect counts for all orbit_indices
        if orbit_indices is None:
            orbit_indices = range(len(self))
        structure_icet = Structure.from_atoms(structure)
        cluster_counts_full = OrderedDict()
        for i in orbit_indices:
            orbit = supercell_orbit_list.get_orbit(i)
            counts = orbit.count_clusters(structure_icet)
            sorted_counts = {}
            for symbols, count in counts.items():
                sorted_counts[tuple(sorted(symbols))] = sorted_counts.get(symbols, 0) + count
            cluster_counts_full[i] = sorted_counts
        return cluster_counts_full
