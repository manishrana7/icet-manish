#!/usr/bin/env python

import unittest

from ase.build import bulk
from icet.core.cluster_counts import ClusterCounts
from icet.core.cluster import Cluster
from icet.core.lattice_site import LatticeSite
from icet.core.neighbor_list import get_neighbor_lists
from icet.core.orbit_list import OrbitList
from icet.core.many_body_neighbor_list import ManyBodyNeighborList
from icet import Structure


def strip_surrounding_spaces(input_string):
    '''
    Helper function that removes both leading and trailing spaces from a
    multi-line string.

    Returns
    -------
    str
        original string minus surrounding spaces and empty lines
    '''
    from io import StringIO
    s = []
    for line in StringIO(input_string):
        if len(line.strip()) == 0:
            continue
        s += [line.strip()]
    return '\n'.join(s)


class TestClusterCounts(unittest.TestCase):
    """
    Container for test of the module functionality.
    """

    def __init__(self, *args, **kwargs):
        super(TestClusterCounts, self).__init__(*args, **kwargs)

        self.atoms = bulk('Ni', 'hcp', a=2.0).repeat([2, 1, 1])
        self.atoms.set_chemical_symbols('NiFeNi2')
        self.cutoffs = [2.2]
        self.neighbor_lists = get_neighbor_lists(
            atoms=self.atoms, cutoffs=self.cutoffs)
        # orbit list
        self.structure = Structure.from_atoms(self.atoms)
        self.orbit_list = OrbitList(self.neighbor_lists,
                                    self.structure)
        self.orbit_list.sort()

    def setUp(self):
        """
        SetUp an empty cluster counts .
        """
        self.cluster_counts = ClusterCounts()

    def test_count_lattice_neighbors(self):
        """
        Test singlet and pair counts given
        many-body neighbor list.
        """
        mbnl = ManyBodyNeighborList()
        mbnl_pairs = mbnl.build(self.neighbor_lists, 0, False)
        self.cluster_counts.count_lattice_neighbors(self.structure, mbnl_pairs)
        cluster_map = self.cluster_counts.get_cluster_counts()

        cluster_singlet = Cluster(self.structure,
                                  [LatticeSite(0, [0., 0., 0.])])
        cluster_pair = Cluster(self.structure,
                               [LatticeSite(0, [0., 0., 0.]),
                                LatticeSite(1, [0., 0., 0.])])
        clusters = [cluster_singlet, cluster_pair]

        expected_counts = [{(28,): 1},
                           {(26, 28): 4, (28, 28): 7}]
        for k, cluster in enumerate(clusters):
            count = cluster_map[cluster]
            self.assertEqual(count, expected_counts[k])

    def test_count_lattice_sites(self):
        """
        Test cluster_counts counts a pair given a set of
        lattice neighbors.
        """
        lattice_sites = []
        lattice_sites.append(LatticeSite(0, [0., 0., 0.]))
        lattice_sites.append(LatticeSite(1, [0., 0., 0.]))

        cluster = Cluster(self.structure, lattice_sites)

        self.cluster_counts.count(self.structure, lattice_sites)
        cluster_map = self.cluster_counts.get_cluster_counts()

        count = cluster_map[cluster]
        self.assertEqual(count, {(26, 28): 1})

    def test_count_list_lattice_sites(self):
        """
        Test cluster_counts counts the right number of pairs
        given a list of lattice neighbors.
        """
        lattice_sites = []
        lattice_sites.append(LatticeSite(0, [0., 0., 0.]))
        lattice_sites.append(LatticeSite(1, [0., 0., 0.]))

        lattice_sites2 = []
        lattice_sites2.append(LatticeSite(0, [0., 0., 0.]))
        lattice_sites2.append(LatticeSite(2, [0., 0., 0.]))

        cluster = Cluster(self.structure, lattice_sites)

        lattice_neighbors = [lattice_sites, lattice_sites2]

        self.cluster_counts.count(self.structure, lattice_neighbors,
                                  cluster, True)
        cluster_map = self.cluster_counts.get_cluster_counts()

        count = cluster_map[cluster]
        self.assertEqual(count, {(28, 26): 1, (28, 28): 1})

    def test_count_orbitlist(self):
        """
        Test cluster_counts given orbits in an orbitlist.
        """
        cluster_singlet = Cluster(self.structure, [], False, 0)
        cluster_pair = Cluster(self.structure, [], False, 1)
        clusters = [cluster_singlet, cluster_pair]

        expected_counts = [{(26,): 1, (28,): 3},
                           {(26, 26): 1, (26, 28): 10, (28, 28): 13}]
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()

        for k, cluster in enumerate(clusters):
            count = cluster_map[cluster]
            self.assertEqual(count, expected_counts[k])

    def test_count_orbitlist_non_pbc(self):
        """
        Test cluster counts using orbitlist for a non-pbc structure.
        """
        atoms_non_pbc = self.atoms.copy()
        atoms_non_pbc.set_pbc(False)
        structure = Structure.from_atoms(atoms_non_pbc)
        neighbor_lists = get_neighbor_lists(atoms=atoms_non_pbc,
                                            cutoffs=self.cutoffs)
        orbit_list = OrbitList(neighbor_lists, structure)

        cluster_singlet = Cluster(self.structure, [], False, 0)
        cluster_pair = Cluster(self.structure, [], False, 1)
        clusters = [cluster_singlet, cluster_pair]

        expected_counts = [{(26,): 1, (28,): 3},
                           {(26, 28): 2, (28, 28): 2}]

        self.cluster_counts.count_clusters(structure,
                                           orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()

        for k, cluster in enumerate(clusters):
            count = cluster_map[cluster]
            self.assertEqual(count, expected_counts[k])

    def test_len(self):
        """
        Test total size of counts.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        self.assertEqual(len(self.cluster_counts),
                         len(self.orbit_list))

    def test_reset(self):
        """
        Test reset.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        self.cluster_counts.reset()
        self.assertEqual(len(self.cluster_counts), 0)

    def test_get_cluster_counts_info(self):
        """
        Test get_cluster_counts_info functionality.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        self.cluster_counts.setup_cluster_counts_info()
        elements, count = self.cluster_counts.get_cluster_counts_info(0)
        self.assertEqual(elements, ['Fe'])
        self.assertEqual(count, 1)

    def test_repr(self):
        """
        Test representation of cluster_counts.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        retval = self.cluster_counts.repr()
        target = '''
Cluster Counts
------------------------------
[0] [] 0.0000
------------------------------
Fe    1
Ni    3
------------------------------
[0, 0] [2.0] 1.0000
------------------------------
Fe   Fe    1
Fe   Ni    10
Ni   Ni    13
'''
        self.assertEqual(strip_surrounding_spaces(target),
                         strip_surrounding_spaces(retval))

    def test_print_overview(self):
        """
        Test print overview.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        self.cluster_counts.print_overview()


if __name__ == '__main__':
    unittest.main()
