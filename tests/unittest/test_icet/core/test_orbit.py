import unittest

from icet.core.lattice_site import LatticeSite
from icet.core.structure import Structure
from icet.core.cluster import Cluster
from icet.core.orbit import Orbit
from ase.build import bulk
import itertools
import numpy as np


class TestOrbit(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestOrbit, self).__init__(*args, **kwargs)

        self.lattice_sites_pairs = []
        indices = [i for i in range(8)]
        unitcell_offsets = []
        cartesian_product_lists = [[0., 1.], [0., 1.], [0., 1.]]
        for element in itertools.product(*cartesian_product_lists):
            unitcell_offsets.append(list(element))
        self.lattice_sites_pairs = [[LatticeSite(index,
                                                 unitcell_offset),
                                     LatticeSite(index + 1,
                                                 unitcell_offset)]
                                    for index, unitcell_offset in
                                    zip(indices, unitcell_offsets)]

        self.lattice_sites_triplets = [[LatticeSite(index,
                                                    unitcell_offset),
                                        LatticeSite(
                                            index + 1, unitcell_offset),
                                        LatticeSite(
                                            index + 3, unitcell_offset)]
                                       for index, unitcell_offset in
                                       zip(indices, unitcell_offsets)]

        self.allowed_permutations_pair = set([(0, 1)])
        self.allowed_permutations_triplet = set([(0, 2, 1)])

    def shortDescription(self):
        """Silences unittest from printing the docstrings in test cases."""
        return None

    def setUp(self):
        """Instantiates class before each test."""
        self.structure = Structure.from_atoms(bulk('Al'))
        lattice_site_for_cluster = [LatticeSite(0, [i, 0, 0]) for i in range(3)]

        self.pair_sites = [lattice_site_for_cluster[0],
                           lattice_site_for_cluster[1]]
        self.triplet_sites = lattice_site_for_cluster

        self.pair_cluster = Cluster(self.structure, self.pair_sites)
        self.triplet_cluster = Cluster(self.structure, self.triplet_sites)

        self.orbit_pair = Orbit(self.structure, [self.pair_sites],
                                self.allowed_permutations_pair)
        self.orbit_triplet = Orbit(self.structure, [self.triplet_sites],
                                   self.allowed_permutations_triplet)

    def test_init(self):
        """Tests the initializer."""
        orbit = Orbit(self.structure, [self.pair_sites], self.allowed_permutations_pair)
        self.assertIsInstance(orbit, Orbit)

        orbit = Orbit(self.structure, [self.triplet_sites], self.allowed_permutations_triplet)
        self.assertIsInstance(orbit, Orbit)

    def test_equivalent_clusters(self):
        """Tests getting the equivalent sites of orbit."""
        self.orbit_pair.equivalent_clusters = self.lattice_sites_pairs
        self.assertEqual(self.orbit_pair.equivalent_clusters,
                         self.lattice_sites_pairs)

    def test_representative_cluster(self):
        """Tests getting the representative cluster of orbit."""
        cluster = self.orbit_pair.representative_cluster
        self.assertEqual(str(cluster), str(self.pair_cluster))

        cluster = self.orbit_triplet.representative_cluster
        self.assertEqual(str(cluster), str(self.triplet_cluster))

    def test_sites_of_representative_cluster(self):
        """Tests getting the representative sites of orbit."""
        self.orbit_pair.equivalent_clusters = self.lattice_sites_pairs
        self.assertEqual(self.orbit_pair.sites_of_representative_cluster,
                         self.lattice_sites_pairs[0])

    def test_order(self):
        """Tests getting the order of orbit."""
        self.assertEqual(
            self.orbit_pair.order, 2)
        self.assertEqual(
            self.orbit_triplet.order, 3)

    def test_len(self):
        """Tests lenght of orbit."""
        self.assertEqual(len(self.orbit_pair), 1)
        self.orbit_pair.equivalent_clusters = self.lattice_sites_pairs
        self.assertEqual(len(self.orbit_pair),
                         len(self.lattice_sites_pairs))

    def test_radius(self):
        """Tests geometrical size of orbit."""
        self.orbit_pair.radius

    def test_add(self):
        """Tests that offset is effectively added to orbit."""
        added_offset = np.array((1., 1., 1.))
        self.orbit_pair.equivalent_clusters = self.lattice_sites_pairs

        # Create a new orbit with offseted equivalent lattice sites
        orbit = self.orbit_pair + added_offset
        self.assertNotEqual(orbit, self.orbit_pair)

        # Loop through sites check index is same and offset has been added
        for sites, sites_with_offset in zip(
                self.orbit_pair.equivalent_clusters,
                orbit.equivalent_clusters):
            for site, site_with_offset in zip(sites, sites_with_offset):
                self.assertEqual(site.index, site_with_offset.index)
                self.assertListEqual(list(site.unitcell_offset),
                                     list(site_with_offset.unitcell_offset -
                                          added_offset))

        # List is allowed but only lenght 3
        with self.assertRaises(TypeError):
            orbit + np.array([1, 1, 1, 1])

    def test_allowed_permutations(self):
        """Tests the allowed permutations property."""
        self.assertEqual(self.orbit_pair.allowed_permutations,
                         [list(i) for i in self.allowed_permutations_pair])
        self.assertEqual(self.orbit_triplet.allowed_permutations,
                         [list(i) for i in self.allowed_permutations_triplet])

    def test_get_mc_vectors_pairs(self):
        """Tests the get mc vectors functionality for a pair orbit."""
        # Binary mc vectors
        # Allow only identity permutation
        allowed_permutations = set([tuple(i for i in range(self.orbit_pair.order))])
        orbit = Orbit(self.structure, [self.pair_sites], allowed_permutations)
        mc_vectors = orbit.get_mc_vectors([2] * orbit.order)
        self.assertEqual(mc_vectors, [[0, 0]])

        # Ternary mc vectors
        mc_vectors = orbit.get_mc_vectors([3] * self.orbit_pair.order)
        target = [[0, 0], [0, 1], [1, 0], [1, 1]]
        self.assertEqual(mc_vectors, target)

        # Allow the [1,0] permutation
        allowed_permutations = set([(0, 1), (1, 0)])
        orbit = Orbit(self.structure, [self.pair_sites], allowed_permutations)
        mc_vectors = orbit.get_mc_vectors([3] * self.orbit_pair.order)
        target = [[0, 0], [0, 1], [1, 1]]
        self.assertEqual(mc_vectors, target)

    def test_get_mc_vectors_triplets(self):
        """Tests  the get mc vectors functionality for a triplet orbit."""
        # Binary mc vectors
        # Allow only identity permutation
        allowed_permutations = set([tuple(i for i in range(self.orbit_triplet.order))])
        orbit = Orbit(self.structure, [self.triplet_sites], allowed_permutations)

        # Ternary mc vectors
        mc_vectors = orbit.get_mc_vectors([3] * self.orbit_triplet.order)
        target = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                  [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
        self.assertEqual(mc_vectors, target)

        # Allow the [0, 2, 1] permutation
        allowed_permutations = set([(0, 1, 2), (0, 2, 1)])
        orbit = Orbit(self.structure, [self.triplet_sites], allowed_permutations)
        mc_vectors = self.orbit_triplet.get_mc_vectors([3] * self.orbit_triplet.order)
        target = [[0, 0, 0], [0, 0, 1], [0, 1, 1],
                  [1, 0, 0], [1, 0, 1], [1, 1, 1]]
        self.assertEqual(mc_vectors, target)


if __name__ == '__main__':
    unittest.main()
