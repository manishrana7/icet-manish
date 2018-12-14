from ase import Atoms
from icet import ClusterExpansion
from mchammer.observers.base_observer import BaseObserver
from icet.core.structure import Structure
from icet.core.local_orbit_list_generator import LocalOrbitListGenerator
from _icet import ClusterCounts as _ClusterCounts
from typing import List
from collections import namedtuple, OrderedDict
from icet.tools.geometry import chemical_number_to_chemical_symbol, chemical_symbols_to_numbers
import pandas as pd
ClusterCountInfo = namedtuple("ClusterCountInfo", ['counts', 'dc_tags'])


class ClusterCountObserver(BaseObserver):
    """
    This class represents a cluster expansion (CE) observer.

    A CE observer allows to compute a property described by a CE along the
    trajectory sampled by a Monte Carlo (MC) simulation. In general this CE
    differs from the CE that is used to generate the trajectory. For example in
    a canonical MC simulation the latter would usually represent an energy
    (total or mixing energy) whereas the former CE(s) could map lattice
    constant or band gap.

    Parameters
    ----------
    cluster_expansion : :class:`icet.ClusterExpansion` cluster expansion model
        to be used for observation
    tag : str
        human readable observer name (default: `ClusterExpansionObserver`)
    interval : int
        observation interval during the Monte Carlo simulation

    Attributes
    ----------
    tag : str
        name of observer
    interval : int
        observation interval
    """

    def __init__(self, cluster_space, atoms: Atoms,
                 interval: int,
                 tag: str = 'ClusterCountObserver') -> None:

        structure = Structure.from_atoms(atoms)
        self._cluster_space = cluster_space
        local_orbit_list_generator = LocalOrbitListGenerator(
            cluster_space.orbit_list, structure)

        self._full_orbit_list = local_orbit_list_generator.generate_full_orbit_list()
        self._cluster_counts_cpp = _ClusterCounts()

        self._cluster_keys = []
        for i, orbit in enumerate(self._full_orbit_list.orbits):
            cluster = orbit.representative_cluster
            cluster.tag = i
            self._cluster_keys.append(cluster)

        super().__init__(interval=interval, return_type=dict, tag=tag)
        self._get_empty_counts()

    def _get_empty_counts(self):
        """Returns the object which will be filled with counts"""
        counts = {}
        for i, cluster in enumerate(self._cluster_keys):
            order = len(cluster)
            dc_tags = []
            possible_decorations = self._cluster_space.get_possible_orbit_decorations(
                cluster.tag)
            for decoration in possible_decorations:
                tag = "{}_{}".format(i, '_'.join(decoration))
                dc_tags.append(tag)
            assert order == len(
                possible_decorations[0]), f"{order} is not {len(possible_decorations[0])}, {possible_decorations}, {cluster}"
            counts_for_this_cluster = {
                decoration: 0 for decoration in possible_decorations}
            count_info = ClusterCountInfo(
                counts=counts_for_this_cluster, dc_tags=dc_tags)

            counts[cluster] = count_info
        return counts

    def _generate_counts(self, atoms):
        """Generates the counts into a pandas dataframe

        Parameters
        ----------
        atoms
            input atomic structure.
        """
        structure = Structure.from_atoms(atoms)
        self._cluster_counts_cpp.count_orbit_list(structure,
                                                  self._full_orbit_list, False, True)

        # self._cluster_counts_cpp.setup_cluster_counts_info()

        empty_counts = self._get_empty_counts()
        observable_dict = {}
        pandas_rows = []

        # std::unordered_map<Cluster, std::map<std::vector<int>, int>>
        cluster_counts = self._cluster_counts_cpp.get_cluster_counts()
        for cluster_key, chemical_number_counts_dict in cluster_counts.items():

            # for chemical_numbers, count in chemical_number_counts_dict.items():
            for chemical_symbols in empty_counts[cluster_key].counts.keys():
                chemical_numbers = tuple(
                    chemical_symbols_to_numbers(chemical_symbols))
                count = chemical_number_counts_dict.get(chemical_numbers, 0)
                pandas_row = {}
                pandas_row['dc_tag'] = "{}_{}".format(
                    cluster_key.tag, '_'.join(chemical_symbols))
                pandas_row['decoration'] = chemical_symbols
                pandas_row['count'] = count
                pandas_row['orbit_index'] = cluster_key.tag
                pandas_row['order'] = len(cluster_key)
                pandas_row['radius'] = cluster_key.radius
                pandas_rows.append(pandas_row)
        self._cluster_counts_cpp.reset()
        self.count_frame = pd.DataFrame(pandas_rows)

    def get_observable(self, atoms: Atoms) -> float:
        """
        Returns the value of the property from a cluster expansion model
        for a given atomic configuration.

        Parameters
        ----------
        atoms
            input atomic structure.
        """
        self._generate_counts(atoms)

        count_dict = {row['dc_tag']: row['count']
                      for i, row in self.count_frame.iterrows()}
        return count_dict
