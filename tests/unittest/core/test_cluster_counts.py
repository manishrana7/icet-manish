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

        self.cutoffs = [1.4, 0.5, 0.8]
        self.atoms_prim = bulk('Ni', 'hcp', a=0.625, c=1.0)
        self.neighbor_list_prim = get_neighbor_lists(
            atoms=self.atoms_prim, cutoffs=self.cutoffs)
        self.structure_prim = Structure.from_atoms(self.atoms_prim)
        self.orbit_list = OrbitList(self.neighbor_list_prim,
                                    self.structure_prim)
        self.orbit_list.sort()
        self.order = []
        self.size = []
        self.multiplicity = []
        for k in range(len(self.orbit_list)):
            orbit = self.orbit_list.get_orbit(k)
            cluster_repr = orbit.get_representative_cluster()
            self.order.append(cluster_repr.order)
            self.size.append(cluster_repr.geometrical_size)
            self.multiplicity.append(len(orbit))

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
        self.assertEqual(len(self.cluster_map), len(self.orbit_list))

    def test_clusters_order(self):
        '''
        Test clusters orders
        '''
        for i, cluster in enumerate(sorted(self.cluster_map.keys())):
            self.assertEqual(cluster.order, self.order[i])

    def test_clusters_size(self):
        '''
        Test clusters size
        '''
        for i, cluster in enumerate(sorted(self.cluster_map.keys())):
            self.assertEqual(cluster.geometrical_size, self.size[i])

    def test_clusters_count(self):
        '''
        Test clusters counts
        '''
        for i, key in enumerate(sorted(self.cluster_map.keys())):
            val = self.cluster_map[key]
            for count in val.values():
                self.assertEqual(count, self.multiplicity[i])

    def test_clusters_count_low_symmetry(self):
        '''
        Test cluster counts is larger for a
        structure with lower symmetry
        TODO: Maybe not an optimal test
        '''
        atoms = self.atoms_prim.copy()
        atoms.set_chemical_symbols(['Ni', 'Fe'])
        structure = Structure.from_atoms(atoms)
        self.cluster_counts.count_clusters(structure, self.orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()
        counts_low_symmetry = [list(val.values())
                               for val in cluster_map.values()]
        counts = [list(val.values()) for val in self.cluster_map.values()]
        self.assertNotEqual(counts, counts_low_symmetry)


if __name__ == '__main__':
    unittest.main()
