#!/usr/bin/env Python3


import unittest

from icet.core_py.orbit import Orbit
from icet.core_py.lattice_site import LatticeSite
from icet.core.lattice_site import LatticeSite as LatticeSite_cpp
from icet.core.cluster import Cluster
from icet.core.orbit import Orbit as Orbit_cpp
from ase.build import bulk
import itertools
import numpy as np


class TestOrbit(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        super(TestOrbit, self).__init__(*args, **kwargs)

        self.lattice_sites_pairs = []
        indices = [i for i in range(8)]
        unitcell_offsets = []
        cartesian_product_lists = [[0., 1.], [0., 1.], [0., 1.]]
        for element in itertools.product(*cartesian_product_lists):
            unitcell_offsets.append(list(element))
        self.lattice_sites_pairs = [[LatticeSite(index, unitcell_offset),
                                     LatticeSite(index + 1, unitcell_offset)]
                                    for index, unitcell_offset in
                                    zip(indices, unitcell_offsets)]
        self.lattice_sites_triplets = [[LatticeSite(index, unitcell_offset),
                                        LatticeSite(
                                            index + 1, unitcell_offset),
                                        LatticeSite(
                                            index + 3, unitcell_offset)]
                                       for index, unitcell_offset in
                                       zip(indices, unitcell_offsets)]
        self.lattice_sites_pairs_cpp = [[LatticeSite_cpp(index,
                                                         unitcell_offset),
                                         LatticeSite_cpp(index + 1,
                                                         unitcell_offset)]
                                        for index, unitcell_offset in
                                        zip(indices, unitcell_offsets)]

        self.lattice_sites_triplets_cpp = [[LatticeSite_cpp(index,
                                                            unitcell_offset),
                                            LatticeSite_cpp(
                                            index + 1, unitcell_offset),
                                            LatticeSite_cpp(
                                            index + 3, unitcell_offset)]
                                           for index, unitcell_offset in
                                           zip(indices, unitcell_offsets)]

    def setUp(self):
        '''
        Instantiate class before each test.
        '''
        self.orbit = Orbit()
        atoms = bulk("Al")
        lattice_site_for_cluster = [
            LatticeSite(0, [i, 0, 0]) for i in range(3)]

        self.pair_cluster = Cluster.from_python(
            atoms, [lattice_site_for_cluster[0],
                    lattice_site_for_cluster[1]], True)
        self.triplet_cluster = Cluster.from_python(
            atoms, lattice_site_for_cluster, True)

        self.orbit_pair_cpp = Orbit_cpp(self.pair_cluster)
        self.orbit_triplet_cpp = Orbit_cpp(self.triplet_cluster)

    def test_init(self):
        '''
        Just testing that the setup
        (initialization) of tested class work
        '''
        # initialize from ASE Atoms
        orbit = Orbit()
        self.assertIsInstance(orbit, Orbit)

        orbit = Orbit_cpp(self.pair_cluster)
        self.assertIsInstance(orbit, Orbit_cpp)

        orbit = Orbit_cpp(self.triplet_cluster)
        self.assertIsInstance(orbit, Orbit_cpp)

    def test_property_equivalent_sites(self):
        '''
        Testing get_orbit_list_info functionality
        '''

        # test Python version
        self.assertEqual(self.orbit.equivalent_sites, [])
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        self.assertEqual(self.orbit.equivalent_sites, self.lattice_sites_pairs)

        # Test C++ version
        self.orbit_pair_cpp.add_equivalent_sites(self.lattice_sites_pairs_cpp)
        self.assertEqual(self.orbit_pair_cpp.get_equivalent_sites(),
                         self.lattice_sites_pairs_cpp)

        # Test C++ and Python together
        self.assertEqual(
            self.orbit_pair_cpp.get_equivalent_sites(),
            self.orbit.equivalent_sites)

    def test_property_representative_cluster(self):
        """
        Test getting the representative Cluster
        TODO
        ----
        * implement this more thoroughly when cluster is added
        """
        self.orbit.representative_cluster
        self.orbit_pair_cpp.get_representative_cluster()
        self.orbit_triplet_cpp.get_representative_cluster()

    def test_property_representative_sites(self):
        """
        Test getting the representative sites.
        """
        with self.assertRaises(IndexError):
            self.orbit.representative_sites
            self.orbit_pair_cpp.get_representative_sites()

        # Test Python version
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        self.assertEqual(self.orbit.representative_sites,
                         self.lattice_sites_pairs[0])

        # Test C++ version
        self.orbit_pair_cpp.add_equivalent_sites(self.lattice_sites_pairs_cpp)
        self.assertEqual(self.orbit_pair_cpp.get_representative_sites(),
                         self.lattice_sites_pairs_cpp[0])

        # Test C++ and Python together
        self.assertEqual(self.orbit_pair_cpp.get_representative_sites(),
                         self.orbit.representative_sites)

    def test_property_order(self):
        """
        Test getting the order from an orbit.
        """
        with self.assertRaises(IndexError):
            self.orbit.order

        self.orbit.equivalent_sites = self.lattice_sites_pairs
        self.assertEqual(self.orbit.order,
                         len(self.lattice_sites_pairs[0]))

    def test_len(self):
        """
        Test len of orbit.
        """
        self.assertEqual(len(self.orbit), 0)
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        self.assertEqual(len(self.orbit), len(self.lattice_sites_pairs))

    def test_property_geometrical_size(self):
        """
        Test the property geometrical size.
        """
        self.orbit.geometrical_size

    def test_sort(self):
        """
        Test the sort functionality.
        """
        self.orbit.equivalent_sites = sorted(self.lattice_sites_pairs,
                                             reverse=True)
        self.orbit.sort()
        self.assertEqual(self.orbit.equivalent_sites,
                         sorted(self.lattice_sites_pairs))

    def test_eq(self):
        """
        Test equality functionality.
        """
        self.assertEqual(self.orbit, self.orbit)

    def test_lt(self):
        """
        Test less than functionality.
        """
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        self.assertFalse(self.orbit < self.orbit)

    def test_add(self):
        """
        Test add operator.
        """
        added_offset = np.array((1., 1., 1.))
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        orbit = self.orbit + added_offset
        self.assertNotEqual(orbit, self.orbit)
        for sites, sites_with_offset in zip(self.orbit.equivalent_sites,
                                            orbit.equivalent_sites):
            for site, site_with_offset in zip(sites, sites_with_offset):
                self.assertEqual(site.index, site_with_offset.index)
                self.assertListEqual(list(site.unitcell_offset),
                                     list(site_with_offset.unitcell_offset -
                                          added_offset))

        with self.assertRaises(TypeError):
            orbit + [1, 1, 1]
        with self.assertRaises(TypeError):
            orbit + np.array([1, 1, 1, 1])

    def test_property_permutations_to_representative(self):
        """
        Test the permutations to representative property.
        """
        self.assertEqual(self.orbit.permutations_to_representative, [])
        self.orbit.permutations_to_representative = [[1, 2, 3]]
        self.assertListEqual(
            self.orbit.permutations_to_representative, [[1, 2, 3]])

    def test_property_allowed_permutations(self):
        """
        Test the allowed permutations property.
        """
        self.assertEqual(self.orbit.allowed_permutations, [])
        self.orbit.allowed_permutations = [[1, 2, 3]]
        self.assertListEqual(self.orbit.allowed_permutations, [[1, 2, 3]])

    def test_property_permutated_sites(self):
        """
        Test the permutated sites property.
        """
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        # Raises IndexError when permutations to primitive is not set
        with self.assertRaises(IndexError):
            self.orbit.permutated_sites

        # Provide the identity permutation
        self.orbit.permutations_to_representative = [
            [i for i in range(self.orbit.order)]] * len(self.orbit)

        self.assertEqual(self.orbit.permutated_sites,
                         self.orbit.equivalent_sites)

        # Provide a completely reversed permutation, [i,j,k] ->[k,j,i]
        self.orbit.permutations_to_representative = [
            [i for i in reversed(range(self.orbit.order))]] * len(self.orbit)

        for perm_sites, sites in zip(self.orbit.permutated_sites,
                                     self.orbit.equivalent_sites):
            self.assertEqual(perm_sites, list(reversed(sites)))

    def test_get_mc_vectors_pairs(self):
        """
        Test  the get mc vectors functionality for a pair orbit
        """
        self.orbit.equivalent_sites = self.lattice_sites_pairs
        # Binary mc vectors
        # Allow only identity permutation
        self.orbit.allowed_permutations = [
            [i for i in range(self.orbit.order)]]
        mc_vectors = self.orbit.get_mc_vectors([2] * self.orbit.order)
        self.assertEqual(mc_vectors, [(0, 0)])

        # Ternary mc vectors
        mc_vectors = self.orbit.get_mc_vectors([3] * self.orbit.order)
        target = [(0, 0), (0, 1), (1, 0), (1, 1)]
        self.assertEqual(mc_vectors, target)

        # Allow the permutation [1,0] permutation
        self.orbit.allowed_permutations = ([0, 1], [1, 0])
        mc_vectors = self.orbit.get_mc_vectors([3] * self.orbit.order)
        target = [(0, 0), (0, 1), (1, 1)]
        self.assertEqual(mc_vectors, target)

    def test_get_mc_vectors_triplets(self):
        """
        Test  the get mc vectors functionality for a triplet orbit
        """
        self.orbit.equivalent_sites = self.lattice_sites_triplets
        # Binary mc vectors
        # Allow only identity permutation
        self.orbit.allowed_permutations = [
            [i for i in range(self.orbit.order)]]

        # Ternary mc vectors
        mc_vectors = self.orbit.get_mc_vectors([3] * self.orbit.order)
        target = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                  (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
        self.assertEqual(mc_vectors, target)

        # Allow the permutation [0,2,1] permutation
        self.orbit.allowed_permutations = ([0, 1, 2], [0, 2, 1])
        mc_vectors = self.orbit.get_mc_vectors([3] * self.orbit.order)
        target = [(0, 0, 0), (0, 0, 1), (0, 1, 1),
                  (1, 0, 0), (1, 0, 1), (1, 1, 1)]
        self.assertEqual(mc_vectors, target)


if __name__ == '__main__':
    unittest.main()
