import unittest
from ase.build import bulk
from icet.core.lattice_site import LatticeSite
from icet.core.cluster import Cluster
from icet.core.orbit import Orbit
from icet.core.orbit_list import OrbitList, create_orbit_list
from icet.core.orbit_list import (
    __get_lattice_site_permutation_matrix as
    get_lattice_site_permutation_matrix)
from icet.core.neighbor_list import get_neighbor_lists
from icet.core.permutation_map import permutation_matrix_from_atoms
from icet.core.structure import Structure
from icet.tools.geometry import get_permutation


class TestOrbitList(unittest.TestCase):
    """Container for test of the module functionality."""
    def __init__(self, *args, **kwargs):
        super(TestOrbitList, self).__init__(*args, **kwargs)
        self.cutoffs = [4.2]
        self.atoms = bulk('Ag', 'sc', a=4.09)

        # representative clusters for testing
        structure = Structure.from_atoms(self.atoms)
        # for singlet
        self.cluster_singlet = Cluster(
            structure, [LatticeSite(0, [0, 0, 0])])
        # for pair
        lattice_sites = [LatticeSite(0, [i, 0, 0]) for i in range(3)]
        self.cluster_pair = Cluster(
            structure, [lattice_sites[0], lattice_sites[1]], True)

    def setUp(self):
        """Instantiate class before each test."""
        # @todo: this can be a single line of code
        permutation_matrix, self.prim_structure, _ = \
            permutation_matrix_from_atoms(self.atoms, self.cutoffs[0])
        self.pm_lattice_sites = \
            get_lattice_site_permutation_matrix(self.prim_structure,
                                                permutation_matrix)
        # @todo: check if the same neighborlist is returned from PermutationMap
        self.neighbor_lists = get_neighbor_lists(
            self.prim_structure, self.cutoffs)

        self.orbit_list = OrbitList(self.prim_structure,
                                    self.pm_lattice_sites,
                                    self.neighbor_lists)

    def test_init(self):
        """Test the different initializers."""
        # empty
        orbit_list = OrbitList()
        self.assertIsInstance(orbit_list, OrbitList)
        # with mnbl and structure
        orbit_list = OrbitList(self.neighbor_lists, self.prim_structure)
        self.assertIsInstance(orbit_list, OrbitList)
        # with mnbl, structure and permutation matrix
        orbit_list = OrbitList(
            self.prim_structure, self.pm_lattice_sites, self.neighbor_lists)
        self.assertIsInstance(self.orbit_list, OrbitList)

    def test_add_orbit(self):
        """Test add orbit funcionality."""
        orbit = Orbit(self.cluster_pair)
        self.orbit_list.add_orbit(orbit)
        self.assertEqual(len(self.orbit_list), 3)

    def test_get_number_of_NClusters(self):
        """Test that only a pair is counted in the orbit list."""
        NPairs = self.orbit_list.get_number_of_NClusters(2)
        self.assertEqual(NPairs, 1)

    def test_get_orbit(self):
        """Test function returns the number of orbits of a given order."""
        # get singlet
        orbit = self.orbit_list.get_orbit(0)
        self.assertEqual(orbit.order, 1)
        # get pair
        orbit = self.orbit_list.get_orbit(1)
        self.assertEqual(orbit.order, 2)
        # check higher order raises an error
        with self.assertRaises(IndexError):
            self.orbit_list.get_orbit(3)

    def test_clear(self):
        """Test orbit list is empty after calling this function."""
        self.orbit_list.clear()
        with self.assertRaises(IndexError):
            self.orbit_list.get_orbit(0)

    def test_sort(self):
        """Test orbits in orbit list are sorted."""
        self.orbit_list.sort()
        for i in range(len(self.orbit_list) - 1):
            self.assertLess(
                self.orbit_list.get_orbit(i), self.orbit_list.get_orbit(i + 1))

    def test_find_orbit(self):
        """
        Test orbit index is retuned from the given representative cluster.
        """
        # @todo: test a non-representative cluster returns -1
        self.assertEqual(
            self.orbit_list.find_orbit(self.cluster_singlet), 0)
        self.assertEqual(
            self.orbit_list.find_orbit(self.cluster_pair), 1)

    def test_is_row_taken(self):
        """Test functionality."""
        taken_rows = set()
        row_indices = tuple([0, 1, 2])
        self.assertFalse(self.orbit_list.is_row_taken(
            taken_rows, row_indices))

        taken_rows = set([row_indices])
        self.assertTrue(self.orbit_list.is_row_taken(
            taken_rows, row_indices))

    def test_get_orbit_list(self):
        """Test a list of orbits is returned from this function."""
        orbit_list = self.orbit_list.get_orbit_list()
        # clusters for testing
        repr_clusters = [self.cluster_singlet, self.cluster_pair]

        for k, orbit in enumerate(orbit_list):
            with self.subTest(orbit=orbit):
                self.assertEqual(orbit.get_representative_cluster(),
                                 repr_clusters[k])

    @unittest.expectedFailure
    def test_get_primitive_structure(self):
        """
        Test get primitive structure functionality.

        Todo
        ----
        Test fails
        """
        self.assertEqual(
            self.orbit_list.get_primitive_structure(), self.prim_structure)

    def test_len(self):
        """Test len of orbit list."""
        self.assertEqual(len(self.orbit_list), 2)

    def test_get_supercell_orbit_list(self):
        """
        Test orbit list is returned for the given supercell

        Todo
        ----
        Test fails for an actual supercell of the testing structure
        """
        atoms_supercell = self.atoms.copy()
        orbit_list_super = \
            self.orbit_list.get_supercell_orbit_list(atoms_supercell)
        orbit_list_super.sort()
        self.orbit_list.sort()
        for k in range(len(orbit_list_super)):
            orbit_super = orbit_list_super.get_orbit(k)
            orbit = self.orbit_list.get_orbit(k)
            self.assertEqual(orbit, orbit_super)

    def test_create_orbit_list(self):
        """
        Test  orbit list is built from structure and cutoffs by calling
        this function.
        """
        orbit_list = create_orbit_list(self.atoms, self.cutoffs)
        for i in range(len(self.orbit_list)):
            orbit = self.orbit_list.get_orbit(i)
            orbit_ = orbit_list.get_orbit(i)
            # check all orbits in both lists are equal
            self.assertEqual(orbit, orbit_)

    def test_equivalent_sites_size(self):
        """
        Test that all the equivalent sites have the same radius.
        """
        atoms = bulk("Al")
        cutoffs = [10, 10]
        orbit_list = create_orbit_list(atoms, cutoffs)
        for orbit in orbit_list.orbits:
            size = orbit.radius
            for eq_sites in orbit.equivalent_sites:
                structure = Structure.from_atoms(atoms)
                cluster = Cluster(structure, eq_sites, True)
                self.assertAlmostEqual(
                    cluster.radius, size, places=5)

    #@unittest.expectedFailure
    def test_allowed_permutations(self):
        """
        Test allowed permutations of orbit.

        Todo
        ----
        Test fails
        """
        atoms = bulk("Al")
        cutoffs = [10, 10]
        orbit_list = create_orbit_list(atoms, cutoffs)
        for orbit in orbit_list.orbits:
            rep_sites = orbit.representative_sites
            translated_sites = orbit_list.get_sites_translated_to_unit_cell(
                rep_sites, False)
            permutations = orbit.allowed_permutations
            for perm in permutations:
                perm_sites = get_permutation(rep_sites, perm)
                self.assertIn(perm_sites, translated_sites)

    def test_orbit_list_non_pbc(self):
        """
        Test that singlets in orbit list retrieves the right number of unique
        sites of the structure with different periodic boundary conditions

        Todo
        ----
        Returned results are incorrect for simple-cubic structures.
        """
        atoms = bulk('Al', 'sc', a=4.0).repeat(4)
        structure = Structure.from_atoms(atoms)
        # [True, True, False]
        structure.set_pbc([True, True, False])
        orbit_list = create_orbit_list(structure, [0.])
        self.assertEqual(len(orbit_list), 4)
        # [True, False, False]
        structure.set_pbc([True, False, False])
        orbit_list = create_orbit_list(structure, [0.])
        self.assertEqual(len(orbit_list), 7)
        # [False]
        structure.set_pbc([False, False, False])
        orbit_list = create_orbit_list(structure, [0.])
        self.assertEqual(len(orbit_list), 20)

    def test_orbit_list_fcc(self):
        """
        Test orbit list has the right number of singlet and pairs for
        a fcc structure.
        """
        atoms = bulk('Al', 'fcc', a=3.0)
        cutoffs = [2.5]
        structure = Structure.from_atoms(atoms)
        orbit_list = create_orbit_list(structure, cutoffs)
        # only a singlet and a pair are expected
        self.assertEqual(len(orbit_list), 2)
        # singlet
        singlet = orbit_list.get_orbit(0)
        self.assertEqual(len(singlet), 1)
        # pair has multiplicity equal to 4
        pairs = orbit_list.get_orbit(1)
        self.assertEqual(len(pairs), 6)
        # not more orbits listed
        with self.assertRaises(IndexError):
            orbit_list.get_orbit(2)

    def test_orbit_list_bcc(self):
        """
        Test orbit list has the right number  of singlet and pairs for
        a bcc structure
        """
        atoms = bulk('Al', 'bcc', a=3.0)
        cutoffs = [3.0]
        structure = Structure.from_atoms(atoms)
        orbit_list = create_orbit_list(structure, cutoffs)
        # one singlet and two pairs expected
        self.assertEqual(len(orbit_list), 3)
        # singlet
        singlet = orbit_list.get_orbit(0)
        self.assertEqual(len(singlet), 1)
        # first pair has multiplicity equal to 4
        pairs = orbit_list.get_orbit(1)
        self.assertEqual(len(pairs), 4)
        # first pair has multiplicity equal to 3
        pairs = orbit_list.get_orbit(2)
        self.assertEqual(len(pairs), 3)
        # not more orbits listed
        with self.assertRaises(IndexError):
            orbit_list.get_orbit(3)

    def test_orbit_list_hcp(self):
        """
        Test orbit list has the right number of singlet and pairs for
        a hcp structure
        """
        atoms = bulk('Ni', 'hcp', a=3.0)
        cutoffs = [3.1]
        structure = Structure.from_atoms(atoms)
        orbit_list = create_orbit_list(structure, cutoffs)
        # only one singlet and one pair expected
        self.assertEqual(len(orbit_list), 3)
        # singlet
        singlet = orbit_list.get_orbit(0)
        self.assertEqual(len(singlet), 2)
        # pair has multiplicity equal to 4
        pairs = orbit_list.get_orbit(1)
        self.assertEqual(len(pairs), 6)
        # pair has multiplicity equal to 4
        pairs = orbit_list.get_orbit(2)
        self.assertEqual(len(pairs), 6)
        # not more orbits listed
        with self.assertRaises(IndexError):
            orbit_list.get_orbit(3)


if __name__ == '__main__':
    unittest.main()
