from collections import namedtuple

import pandas as pd

from _icet import ClusterCounts as _ClusterCounts
from ase import Atoms
from icet.core.local_orbit_list_generator import LocalOrbitListGenerator
from icet.core.structure import Structure
from icet.tools.geometry import chemical_symbols_to_numbers
from mchammer.observers.base_observer import BaseObserver

ClusterCountInfo = namedtuple("ClusterCountInfo", ['counts', 'dc_tags'])


class ClusterCountObserver(BaseObserver):
    """
    This class represents a cluster count observer.

    A cluster count observer allows to count the decorations of
    clusters along the trajectory sampled by a Monte Carlo (MC)
    simulation. For example several canonical MC simulations could
    be executed at different temperatures and the temperature dependence of
    the number of nearest neigbhors of a particular species could
    accessed with this observer.

    Parameters
    ----------
    cluster_space : :class:`icet.ClusterSpace` cluster space to define
        the cluster to be counted
    interval : int
        observation interval during the Monte Carlo simulation

    Attributes
    ----------
    tag : str
        human readable observer name (`ClusterCountObserver`)
    interval : int
        observation interval
    """

    def __init__(self, cluster_space, atoms: Atoms,
                 interval: int) -> None:

        structure = Structure.from_atoms(atoms)
        self._cluster_space = cluster_space
        local_orbit_list_generator = LocalOrbitListGenerator(
            cluster_space.orbit_list, structure)

        self._full_orbit_list = \
            local_orbit_list_generator.generate_full_orbit_list()
        self._cluster_counts_cpp = _ClusterCounts()

        self._cluster_keys = []
        for i, orbit in enumerate(self._full_orbit_list.orbits):
            cluster = orbit.representative_cluster
            cluster.tag = i
            self._cluster_keys.append(cluster)

        super().__init__(interval=interval, return_type=dict,
                         tag='ClusterCountObserver')
        self._get_empty_counts()

    def _get_empty_counts(self):
        """Returns the object which will be filled with counts"""
        counts = {}
        for i, cluster in enumerate(self._cluster_keys):
            order = len(cluster)
            dc_tags = []
            possible_decorations =\
                self._cluster_space.get_possible_orbit_decorations(
                    cluster.tag)
            for decoration in possible_decorations:
                tag = "{}_{}".format(i, '_'.join(decoration))
                dc_tags.append(tag)
            assert order == len(
                possible_decorations[0]), "{} is not {}, {}, {}".format(
                    order, len(possible_decorations[0]), possible_decorations)
            counts_for_this_cluster = {
                decoration: 0 for decoration in possible_decorations}
            count_info = ClusterCountInfo(
                counts=counts_for_this_cluster, dc_tags=dc_tags)

            counts[cluster] = count_info
        return counts

    def _generate_counts(self, atoms: Atoms) ->None:
        """Generates the counts into a pandas dataframe

        Parameters
        ----------
        atoms
            input atomic structure.
        """
        structure = Structure.from_atoms(atoms)
        self._cluster_counts_cpp.count_orbit_list(
            structure, self._full_orbit_list, False, True)

        empty_counts = self._get_empty_counts()
        pandas_rows = []

        # std::unordered_map<Cluster, std::map<std::vector<int>, int>>
        cluster_counts = self._cluster_counts_cpp.get_cluster_counts()
        for cluster_key, chemical_number_counts_dict in \
                cluster_counts.items():

            for chemical_symbols in empty_counts[cluster_key].counts.keys():
                chemical_numbers = tuple(
                    chemical_symbols_to_numbers(chemical_symbols))
                count = chemical_number_counts_dict.get(chemical_numbers, 0)
                pandas_row = {}
                pandas_row['dc_tag'] = "{}_{}".format(
                    cluster_key.tag, '_'.join(chemical_symbols))
                pandas_row['decoration'] = chemical_symbols
                pandas_row['cluster_count'] = count
                pandas_row['orbit_index'] = cluster_key.tag
                pandas_row['order'] = len(cluster_key)
                pandas_row['radius'] = cluster_key.radius
                pandas_rows.append(pandas_row)
        self._cluster_counts_cpp.reset()
        self.count_frame = pd.DataFrame(pandas_rows)

    def get_observable(self, atoms: Atoms) -> dict:
        """
        Returns the value of the property from a cluster expansion model
        for a given atomic configuration.

        Parameters
        ----------
        atoms
            input atomic structure.
        """
        self._generate_counts(atoms)

        count_dict = {row['dc_tag']: row['cluster_count']
                      for i, row in self.count_frame.iterrows()}
        return count_dict
