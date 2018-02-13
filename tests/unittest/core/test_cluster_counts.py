#!/usr/bin/env python

import unittest

from ase.build import bulk
from icet.core.cluster_counts import ClusterCounts
from icet.core.neighbor_list import get_neighbor_lists
from icet.core.orbit_list import OrbitList
from icet.core.many_body_neighbor_list import ManyBodyNeighborList
from icet import Structure


class TestClusterCounts(unittest.TestCase):
    """
    Container for test of the module functionality.

    WIP: Uncompleted tests.
    """

    def __init__(self, *args, **kwargs):
        super(TestClusterCounts, self).__init__(*args, **kwargs)

        self.atoms = bulk('Ni', 'hcp', a=0.625, c=1.0).repeat(2)
        self.neighbor_lists = get_neighbor_lists(
            atoms=self.atoms, cutoffs=[1.4, 0.5, 0.8])
        # orbit list
        self.structure = Structure.from_atoms(self.atoms)
        self.orbit_list = OrbitList(self.neighbor_lists,
                                    self.structure)
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
        """
        SetUp an empty cluster counts .
        """
        self.cluster_counts = ClusterCounts()

    def test_count_lattice_neighbor(self):
        """
        Test lattice neighbors.
        """
        mbnl = ManyBodyNeighborList()
        mbnl_pairs = mbnl.build(self.neighbor_lists, 0, False)
        self.cluster_counts.count_lattice_neighbors(self.structure, mbnl_pairs)
        cluster_map = self.cluster_counts.get_cluster_counts()

        self.assertEqual(len(cluster_map), len(self.orbit_list))

    def test_count(self):
        """
        Test count cluster.
        """
        pass

    def test_count_orbitlist(self):
        """
        Test count cluster.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()

        self.assertEqual(len(cluster_map), len(self.orbit_list))

    def test_clusters_count(self):
        """
        Test multiplicity, order and size of clusters retrieved
        from clusters count.
        """
        self.cluster_counts.count_clusters(self.structure,
                                           self.orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()

        for i, cluster in enumerate(sorted(cluster_map.keys())):
            self.assertEqual(cluster.order, self.order[i])
            self.assertEqual(cluster.geometrical_size, self.size[i])
            counts = cluster_map[cluster]
            for count in counts.values():
                self.assertEqual(count, self.multiplicity[i])

    def test_clusters_count_low_symmetry(self):
        """
        Test cluster counts is larger for a
        structure with lower symmetry
        """
        atoms = self.atoms.copy()
        # atoms.set_chemical_symbols(['Ni', 'Fe'])
        structure = Structure.from_atoms(atoms)
        cluster_counts = ClusterCounts()
        cluster_map_low_symm = cluster_counts.get_cluster_counts()

        cluster_counts.count_clusters(structure, self.orbit_list, False)
        cluster_map = self.cluster_counts.get_cluster_counts()
        counts_low_symm = [list(val.values())
                           for val in cluster_map_low_symm.values()]
        counts = [list(val.values()) for val in cluster_map.values()]
        self.assertEqual(counts, counts_low_symm)


if __name__ == '__main__':
    unittest.main()
