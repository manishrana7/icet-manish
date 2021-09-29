from typing import Dict, List

import pandas as pd

from ase import Atoms
from icet.core.local_orbit_list_generator import LocalOrbitListGenerator
from icet.core.structure import Structure
from mchammer.observers.base_observer import BaseObserver


class ClusterCountObserver(BaseObserver):
    """
    This class represents a cluster count observer.

    A cluster count observer enables one to keep track of the
    occupation of clusters along the trajectory sampled by a Monte
    Carlo (MC) simulation. For example, using this observer, several
    canonical MC simulations could be carried out at different
    temperatures and the temperature dependence of the number of
    nearest neigbhors of a particular species could be accessed with
    this observer.

    Parameters
    ----------
    cluster_space : icet.ClusterSpace
     cluster space to define the clusters to be counted
    structure : ase.Atoms
        defines the lattice that the observer will work on
    interval : int
        observation interval during the Monte Carlo simulation
    max_orbit : int
        only include orbits up to the orbit with this index
        (default is to include all orbits)

    Attributes
    ----------
    tag : str
        human readable observer name
    interval : int
        the observation interval, defaults to None meaning that if the
        observer is used in a Monte Carlo simulation, then the Ensemble object
        will set the interval.
    """

    def __init__(self, cluster_space, structure: Atoms,
                 interval: int = None,
                 orbit_indices: List[int] = None) -> None:
        super().__init__(interval=interval, return_type=dict, tag='ClusterCountObserver')

        self._cluster_space = cluster_space
        local_orbit_list_generator = LocalOrbitListGenerator(
            orbit_list=cluster_space.orbit_list,
            structure=Structure.from_atoms(structure),
            fractional_position_tolerance=cluster_space.fractional_position_tolerance)

        self._full_orbit_list = local_orbit_list_generator.generate_full_orbit_list()

        if orbit_indices is None:
            self._orbit_indices = list(range(len(self._full_orbit_list)))
        elif not isinstance(orbit_indices, List):
            raise ValueError('Argument orbit_indices should be a list of integers, '
                             f'not {type(orbit_indices)}')
        else:
            self._orbit_indices = orbit_indices

        self._empty_counts = self._get_empty_counts()

    def _get_empty_counts(self) -> Dict[int, Dict[List[str], int]]:
        """ Returns the object which will be filled with counts. """
        counts = {}
        for i in self._orbit_indices:
            order = self._full_orbit_list.get_orbit(i).order
            possible_occupations = self._cluster_space.get_possible_orbit_occupations(i)
            assert order == len(possible_occupations[0]), '{} is not {}, {}'.format(
                order, len(possible_occupations[0]), possible_occupations)

            counts[i] = {occupation: 0 for occupation in possible_occupations}
        return counts

    def _generate_counts(self, structure: Atoms) -> None:
        """Counts the occurrence of different clusters and stores this
        information in a pandas dataframe.

        Parameters
        ----------
        structure
            input atomic structure.
        """
        pandas_rows = []
        structure_icet = Structure.from_atoms(structure)
        for i in self._orbit_indices:
            orbit = self._full_orbit_list.get_orbit(i)
            cluster_counts = self._full_orbit_list.get_orbit(i).count_clusters(structure_icet)
            for chemical_symbols in self._empty_counts[i].keys():
                count = cluster_counts.get(chemical_symbols, 0)
                pandas_row = {}
                pandas_row['dc_tag'] = '{}_{}'.format(i, '_'.join(chemical_symbols))
                pandas_row['occupation'] = chemical_symbols
                pandas_row['cluster_count'] = count
                pandas_row['orbit_index'] = i
                pandas_row['order'] = orbit.order
                pandas_row['radius'] = orbit.radius
                pandas_rows.append(pandas_row)
        self.count_frame = pd.DataFrame(pandas_rows)

    def get_observable(self, structure: Atoms) -> dict:
        """
        Returns the value of the property from a cluster expansion model
        for a given atomic configuration.

        Parameters
        ----------
        structure
            input atomic structure
        """
        self._generate_counts(structure)

        count_dict = {row['dc_tag']: row['cluster_count']
                      for i, row in self.count_frame.iterrows()}
        return count_dict
