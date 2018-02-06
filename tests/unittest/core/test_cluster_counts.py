#!/usr/bin/env python

import unittest

from ase.build import bulk
from icet.core.cluster_counts import ClusterCounts
from icet.core.neighbor_list import get_neighbor_lists
from icet.core.orbit_list import OrbitList
from icet import Structure


class TestClusterCounts(unittest.TestCase):
    '''
    Container for test of the module functionality.
    '''

    def __init__(self, *args, **kwargs):
        super(TestClusterCounts, self).__init__(*args, **kwargs)
        self.counts = []
        self.distances = []
        self.order = []

        self.cutoffs = [1.4, 0.5, 0.8]
        self.atoms_prim = bulk('Ni', 'hcp', a=0.625, c=1.0)
        self.neighbor_list_prim = get_neighbor_lists(
            atoms=self.atoms_prim, cutoffs=self.cutoffs)
        self.structure_prim = Structure.from_atoms(self.atoms_prim)
        self.orbit_list = OrbitList(self.neighbor_list_prim,
                                    self.structure_prim)
        self.orbit_list.sort()

    def setUp(self):
        '''
        SetUp.
        '''
        self.cluster_counts = ClusterCounts()
        self.cluster_counts.count_clusters(self.structure_prim,
                                           self.orbit_list, False)
        self.cluster_map = self.cluster_counts.get_cluster_counts()

    def test_len_cluster_map(self):
        '''
        Test size of cluster map
        '''
        pass

    def test_clusters_order(self):
        '''
        Test clusters orders
        '''
        pass

    def test_clusters_distances(self):
        '''
        Test clusters distances
        '''
        pass

    def test_clusters_counts(self):
        '''
        Test clusters counts
        '''
        pass


if __name__ == '__main__':
    unittest.main()
