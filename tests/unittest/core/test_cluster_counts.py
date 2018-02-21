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

    WIP: Uncompleted tests.
    """

    def __init__(self, *args, **kwargs):
        super(TestClusterCounts, self).__init__(*args, **kwargs)

        self.atoms = bulk('Al', 'fcc', cubic=True, a=3.0)
        self.atoms.set_chemical_symbols('AlCuAl2')
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

    def test_count_lattice_neighbor(self):
        """
        Test cluster counts from a many-body neighbor list.
        """
        mbnl = ManyBodyNeighborList()
        mbnl_pairs = mbnl.build(self.neighbor_lists, 0, False)
        self.cluster_counts.count_lattice_neighbors(self.structure, mbnl_pairs)
        cluster_map = self.cluster_counts.get_cluster_counts()

        expected_counts = [{(13,): 1},
                           {(13, 13): 8, (13, 29): 4}]
        for k, count in enumerate(cluster_map.values()):
            self.assertEqual(count, expected_counts[k])

    def test_count_lattice_sites(self):
        """
        Test count clusters thorugh lattice neighbors.
        """
        lattice_sites = []
        lattice_sites.append(LatticeSite(0, [0., 0., 0.]))
        lattice_sites.append(LatticeSite(1, [0., 0., 0.]))

        cluster = Cluster(self.structure, lattice_sites)

        self.cluster_counts.count(self.structure, lattice_sites)
        cluster_map = self.cluster_counts.get_cluster_counts()
        for key, val in cluster_map.items():
            self.assertEqual(key, cluster)
            self.assertTrue(val == {(13, 29): 1})

    def test_count_list_lattice_sites(self):
        """
        Test cluster count thorugh a list of lattice neighbors.
        """
        lattice_sites = []
        lattice_sites.append(LatticeSite(0, [0., 0., 0.]))
        lattice_sites.append(LatticeSite(1, [0., 0., 0.]))

        lattice_sites2 = []
        lattice_sites2.append(LatticeSite(0, [0., 0., 0.]))
        lattice_sites2.append(LatticeSite(2, [0., 0., 0.]))

        cluster = Cluster(self.structure, lattice_sites)

        lattice_neighbors = [lattice_sites, lattice_sites2]

        self.cluster_counts.count(self.structure, lattice_neighbors, cluster)
        cluster_map = self.cluster_counts.get_cluster_counts()
        for key, val in cluster_map.items():
            self.assertEqual(key, cluster)
            self.assertTrue(val == {(13, 13): 1, (13, 29): 1})

    def test_count_orbitlist(self):
        """
        Test cluster count using orbitlist.
        @todo: It might be also good to test clusters but
               clusters retrieved from this count seems to be
               not sorted.
        """
        expected_counts = [{(13,): 3, (29,): 1},
                           {(13, 13): 12, (13, 29): 12}]
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()

        for k, count in enumerate(cluster_map.values()):
            self.assertEqual(count, expected_counts[k])

    def test_count_orbitlist_non_pbc(self):
        """
        Test cluster count using orbotlist for a non-pbc structure.
        """
        atoms_non_pbc = self.atoms.copy()
        atoms_non_pbc.set_pbc(False)
        structure = Structure.from_atoms(atoms_non_pbc)
        neighbor_lists = get_neighbor_lists(atoms=atoms_non_pbc,
                                            cutoffs=self.cutoffs)
        orbit_list = OrbitList(neighbor_lists, structure)

        expected_counts = [{(13,): 3, (29,): 1},
                           {(13, 13): 3, (13, 29): 3}]

        self.cluster_counts.count_clusters(structure,
                                           orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()
        for k, count in enumerate(cluster_map.values()):
            self.assertEqual(count, expected_counts[k])

    def test_reset(self):
        """
        Test reset.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        self.cluster_counts.reset()
        self.assertEqual(len(self.cluster_counts), 0)

    def test_get_cluster_count_info(self):
        """
        Test get cluster count info functionality.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        self.cluster_counts.setup_cluster_counts_info()
        elements, count = self.cluster_counts.get_cluster_counts_info(0)
        self.assertEqual(elements, ['Al'])
        self.assertEqual(count, 3)

    def test_repr(self):
        """
        Test representation.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        retval = self.cluster_counts.repr()
        target = '''
Cluster Counts
------------------------------
[0] [] 0.0000
------------------------------
Al    3
Cu    1
------------------------------
[0, 0] [2.12132] 1.0607
------------------------------
Al   Al    12
Al   Cu    12
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
