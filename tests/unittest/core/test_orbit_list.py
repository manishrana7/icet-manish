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


class TestOrbitList(unittest.TestCase):
    """
    Container for test of the module functionality.
    """

    def __init__(self, *args, **kwargs):
        super(TestOrbitList, self).__init__(*args, **kwargs)
        self.cutoffs = [4.2]
        self.atoms = bulk('Ag', 'sc', a=4.09)

    def setUp(self):
        """ Instantiate class before each test. """
        permutation_matrix, self.prim_structure, _ = \
            permutation_matrix_from_atoms(self.atoms, self.cutoffs[0])
        pm_lattice_sites = \
            get_lattice_site_permutation_matrix(self.prim_structure,
                                                permutation_matrix)
        self.neighbor_lists = get_neighbor_lists(
            self.prim_structure, self.cutoffs)

        self.orbit_list = OrbitList(self.prim_structure,
                                    pm_lattice_sites,
                                    self.neighbor_lists)

    def test_init(self):
        """ Test the different initializers. """
        # empty
        orbit_list = OrbitList()
        self.assertIsInstance(orbit_list, OrbitList)
        # initilize with mbnl
        orbit_list = OrbitList(self.neighbor_lists, self.prim_structure)
        self.assertIsInstance(orbit_list, OrbitList)
        # initiliaze with permutation matrix
        self.assertIsInstance(self.orbit_list, OrbitList)

    def test_add_orbit(self):
        """ Test add orbit funcionality. """
        lattice_site_for_cluster = [
            LatticeSite(0, [i, 0, 0]) for i in range(2)]
        pair_cluster = Cluster.from_python(
            self.atoms, [lattice_site_for_cluster[0],
                         lattice_site_for_cluster[1]], True)
        orbit = Orbit(pair_cluster)
        self.orbit_list.add_orbit(orbit)
        self.assertEqual(len(self.orbit_list), 3)

    def test_get_number_of_NClusters(self):
        """ Test that only a pair is counted in the orbit list. """
        NPairs = self.orbit_list.get_number_of_NClusters(2)
        self.assertEqual(NPairs, 1)

    def test_get_orbit(self):
        """ Test function returns the number of orbits of a given order. """
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
        """ Test orbit list is empty after calling this function. """
        self.orbit_list.clear()
        with self.assertRaises(IndexError):
            self.orbit_list.get_orbit(0)

    def test_sort(self):
        """ Test orbits in orbit list are sorted. """
        self.orbit_list.sort()
        for i in range(len(self.orbit_list) - 1):
            self.assertLess(
                self.orbit_list.get_orbit(i), self.orbit_list.get_orbit(i + 1))

    def test_get_orbit_list(self):
        """ Test getter for orbit list (depends on initializer). """
        lattice_sites = \
            [[LatticeSite(0, [0., 0., 0.])],
             [LatticeSite(0, [0., 0., 0.]), LatticeSite(0, [-1., 0., 0.])],
             [LatticeSite(0, [-1., 1., -1.]), LatticeSite(0, [0., 0., 0.])]]
        orbit_list = self.orbit_list.get_orbit_list()

        for i in range(len(self.orbit_list)):
            self.assertEqual(orbit_list[i].representative_sites,
                             lattice_sites[i])

    def test_get_primitive_structure(self):
        """ Test primitive structure. """
        prim_structure = self.orbit_list.get_primitive_structure()
        self.assertEqual(prim_structure.positions.tolist(),
                         self.prim_structure.positions.tolist())

    def test_len(self):
        """ Test len of orbit list. """
        self.assertEqual(len(self.orbit_list), 2)

    def test_get_supercell_orbit_list(self):
        """
        Test orbit list is returned for the given supercell.

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
        Test  orbit list is built from structure and cutoffs through
        this function.
        """
        structure = Structure.from_atoms(self.atoms)
        orbit_list = create_orbit_list(structure, self.cutoffs)
        for i in range(len(self.orbit_list)):
            orbit = self.orbit_list.get_orbit(i)
            orbit_ = orbit_list.get_orbit(i)
            # all orbits should be equal
            self.assertEqual(orbit, orbit_)

    def test_orbit_list_non_pbc(self):
        """
        Test that singlets in orbit list retrieves the right number of unique
        sites of the structure with different periodic boundary conditions.

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
        a bcc structure.
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
        a hcp structure.
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
