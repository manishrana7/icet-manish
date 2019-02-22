from collections import namedtuple
from typing import Dict

import numpy as np

from ase import Atoms
from icet import ClusterSpace
from mchammer.observers import ClusterCountObserver
from mchammer.observers.base_observer import BaseObserver

ClusterCountInfo = namedtuple("ClusterCountInfo", ['counts', 'dc_tags'])


class BinaryShortRangeOrderObserver(BaseObserver):
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
    cluster_space : :class:`icet.ClusterSpace` cluster space used
    for initialization
    interval : int
        observation interval during the Monte Carlo simulation
    radius : float
        the maximum radius  for the neigbhor shells considered
    Attributes
    ----------
    tag : str
        human readable observer name (`ClusterExpansionObserver`)
    interval : int
        observation interval



    Approach:
    Find all sublattices, A, B, C etc...
    Example:

    A: Pd, Au, ...
    B: H, He, ...

    Map this
    Pd -> A
    Au -> A

    H -> B
    He -> B


    Also know the sites for the different lattices:
    A : [0, 2, 4, 6]
    B : [1, 3, 5, 7]


    Implement
    _get_concentrations() -> Dict[str, float]
    returns concentrations for each element
    """

    def __init__(self, cluster_space, structure: Atoms,
                 interval: int, radius: float) -> None:

        self._structure = structure

        self._cluster_space = ClusterSpace(
            atoms=cluster_space.primitive_structure,
            cutoffs=[radius],
            chemical_symbols=cluster_space.chemical_symbols)
        print(self._cluster_space)
        self._cluster_count_observer = ClusterCountObserver(
            cluster_space=self._cluster_space, atoms=structure,
            interval=interval)

        self._sublattices = self._cluster_space.get_sublattices(structure)
        binary_sublattice_counts = 0
        for symbols in self._sublattices.allowed_species:
            if len(symbols) == 2:
                binary_sublattice_counts += 1
                self._symbols = sorted(symbols)
            elif len(symbols) > 2:
                raise ValueError("Cluster space have more than two allowed"
                                 " species on a sublattice. "
                                 "Allowed species: {}".format(symbols))
        if binary_sublattice_counts != 1:
            raise ValueError("Number of binary sublattices must equal one,"
                             " not {}".format(binary_sublattice_counts))

        super().__init__(interval=interval, return_type=dict,
                         tag='BinaryShortRangeOrderObserver')

    def get_observable(self, atoms: Atoms) -> Dict[str, float]:
        """
        Returns the value of the property from a cluster expansion
        model for a given atomic configuration.

        Parameters
        ----------
        atoms
            input atomic structure.
        """
        self._cluster_count_observer._generate_counts(atoms)
        df = self._cluster_count_observer.count_frame
        print(df['orbit_index'].tolist())
        pair_orbit_indices = set(
            df.loc[df['order'] == 2]['orbit_index'].tolist())

        sro_parameters = {}
        for k, orbit_index in enumerate(sorted(pair_orbit_indices)):
            orbit_df = df.loc[df['orbit_index'] == k]
            A_B_pair_count = 0
            total_count = 0
            for i, row in orbit_df.iterrows():
                total_count += row.cluster_count
                if self._symbols[0] in row.decoration and \
                        self._symbols[1] in row.decoration:
                    A_B_pair_count += row.cluster_count

            key = 'sro_{}_{}'.format(self._symbols[0], k+1)
            conc_B = self._get_concentrations(atoms)[self._symbols[1]]
            value = 1 - A_B_pair_count/(total_count*(1-conc_B))
            sro_parameters[key] = value

        return sro_parameters

    def _get_concentrations(self, structure: Atoms) -> Dict[str, float]:
        """ Returns concentrations for each species in the structure

            Parameters
            ----------
            structure
                the configuration that the concentration will be
                calculated on
        """
        decoration = np.array(structure.get_chemical_symbols())
        concentrations = {}
        for symbols in self._sublattices.allowed_species:
            for symbol in symbols:
                sublattice_index = self._sublattices.get_sublattice_index(
                    symbol=symbol)
                indices = self._sublattices.get_sublattice_sites(
                    sublattice_index)
                symbol_count = decoration[indices].tolist().count(symbol)
                concentration = symbol_count / len(indices)
                concentrations[symbol] = concentration
        return concentrations
