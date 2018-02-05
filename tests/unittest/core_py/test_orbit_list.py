#!/usr/bin/env Python3


import unittest

from icet.core_py.orbit_list import OrbitList
from icet.core_py.lattice_site import LatticeSite
from icet.core_py.permutation_matrix import PermutationMatrix
# from icet.core.orbit_list import create_orbit_list
from ase.build import bulk
# from icet import Structure
from icet import ClusterSpace


class TestOrbitList(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        super(TestOrbitList, self).__init__(*args, **kwargs)
        self.cutoffs = [4.2]
        self.atoms = bulk('Ag', a=4.09)

    def setUp(self):
        '''
        Instantiate class before each test.
        '''
        self.orbit_list = OrbitList(self.atoms, self.cutoffs)
        # self._structure = Structure.from_atoms(self.atoms)

        self.cluster_space_cpp = ClusterSpace(
            self.atoms, self.cutoffs, ["Al", "H"])

    def test_init(self):
        '''
        Just testing that the setup
        (initialization) of tested class work
        '''
        # initialize from ASE Atoms
        orbit_list = OrbitList(self.atoms, self.cutoffs)
        self.assertIsInstance(orbit_list, OrbitList)

    def test_sort(self):
        '''
        Testing len functionality
        '''
        self.orbit_list.sort()
        for i in range(len(self.orbit_list) - 1):
            self.assertLess(
                self.orbit_list.orbits[i], self.orbit_list.orbits[i + 1])

    def test_property_primitive_structure(self):
        '''
        Testing get_orbit_list_info functionality
        '''
        self.orbit_list.primitive_structure
        self.assertEqual(
            self.orbit_list.primitive_structure,
            self.orbit_list.permutation_matrix.primitive_structure)

    def test_property_orbit(self):
        """
        Test orbit property.
        """
        self.orbit_list.orbits
        self.assertEqual(len(self.orbit_list), len(self.orbit_list.orbits))

    def test_is_new_orbit(self):
        """
        Test is new orbit method
        Test this method by going through
        all sites in the orbits and testing
        that there are no new orbits
        that could be made from these sites,
        i.e. they are already taken.
        """
        for orbit in self.orbit_list.orbits:
            for sites in orbit.equivalent_sites:
                self.assertFalse(self.orbit_list.is_new_orbit(sites))
        

    def test_make_orbit(self):
        """
        Test make a new orbit.
        """
        pass

    def test_get_rows(self):
        """
        Test the get row method.
        """
        pass

    def test_get_indices(self):
        """
        Test the get indices method
        """
        pass

    def test_get_all_translated_sites(self):
        """
        Test teh get all translated sites functionality.
        """
        sites = [LatticeSite(0, [0, 0, 0])]
        target = [[LatticeSite(0, [0, 0, 0])]]
        self.assertListEqual(
            self.orbit_list.get_all_translated_sites(sites), target)

        # Does it break when the offset is floats?
        sites = [LatticeSite(0, [0.0, 0.0, 0.0])]
        target = [[LatticeSite(0, [0.0, 0.0, 0.0])]]
        self.assertListEqual(
            self.orbit_list.get_all_translated_sites(sites), target)

        sites = [LatticeSite(0, [1.0, 0.0, 0.0])]
        target = [[LatticeSite(0, [1.0, 0.0, 0.0])], [
            LatticeSite(0, [0.0, 0.0, 0.0])]]
        self.assertListEqual(
            self.orbit_list.get_all_translated_sites(sites), target)

        sites = [LatticeSite(0, [1.0, 0.0, 0.0]),
                 LatticeSite(0, [0.0, 0.0, 0.0])]
        target = [[LatticeSite(0, [1.0, 0.0, 0.0]),
                   LatticeSite(0, [0.0, 0.0, 0.0])],
                  [LatticeSite(0, [0.0, 0.0, 0.0]),
                   LatticeSite(0, [-1, 0.0, 0.0])]]
        self.assertListEqual(
            self.orbit_list.get_all_translated_sites(sites), target)

        sites = [LatticeSite(0, [1.0, 2.0, -1.0]),
                 LatticeSite(2, [2.0, 0.0, 0.0])]

        target = [sites,
                  [LatticeSite(0, [0.0, 0.0, 0.0]),
                   LatticeSite(2, [1.0, -2.0, 1.0])],
                  [LatticeSite(0, [-1.0, 2.0, -1.0]),
                   LatticeSite(2, [0.0, 0.0, 0.0])]]
        self.assertListEqual(
            self.orbit_list.get_all_translated_sites(sites), target)

    def test_no_duplicates_in_orbit(self):
        """
        Test that there are no duplicate sites in orbits
        equivalent sites.
        """
        for orbit in self.orbit_list.orbits:
            for i, site_i in enumerate(orbit.equivalent_sites):
                for j, site_j in enumerate(orbit.equivalent_sites):
                    if j <= i:
                        continue
                    self.assertNotEqual(sorted(site_i), sorted(site_j))

    def test_property_permutation_matrix(self):
        '''
        Test the permutation matrix property.
        '''
        self.assertIsInstance(
            self.orbit_list.permutation_matrix, PermutationMatrix)
        self.assertEqual(
            self.orbit_list.permutation_matrix.cutoff, max(self.cutoffs))

    def test_str(self):
        print(self.orbit_list)
        print(self.cluster_space_cpp)

        for orbit in self.orbit_list.orbits:
            for sites in sorted(orbit.equivalent_sites):
                for site in sorted(sites):
                    print(site, end=' ')
                print()
            print("----")
        print("C++ version")
        for orbit in self.cluster_space_cpp.get_orbit_list().get_orbit_list():
            for sites in sorted(orbit.equivalent_sites):
                for site in sorted(sites):
                    print(site, end=' ')
                print()
            print("----")

    def test_singlets_particle(self):
        """
        Test that a particle get correct number of singlets.
        """
        #  One corner site,X,
        #  one next to the corner along the side, O,
        #  one site on the surface, K,
        #  and one of the center 2x2x2 sites J
        #  and in the darkness to bind them all
        #
        #    from the side:
        #    X O O X
        #    H K K H
        #    H K K H
        #    X O O X
        #
        #  below the K sites are another J site in the 2x2x2 positions
        #
        atoms = bulk("Al", 'sc', a=1).repeat(4)
        atoms.pbc = False

        orbit_list = OrbitList(atoms, [0])
        cluster_space = ClusterSpace(atoms, [0], ["Al", "W"])
        self.assertEqual(len(orbit_list), 4)

        # Size 3 also has 4 unique sites
        atoms = bulk("Al", 'sc', a=1).repeat(3)
        atoms.pbc = False

        orbit_list = OrbitList(atoms, [0])
        cluster_space = ClusterSpace(atoms, [0], ["Al", "W"])
        self.assertEqual(len(orbit_list), 4)

        # Test partial pbc:
        # Now we only distinguish sites by layers from
        # The surface. We have 4 layers but with mirror
        # symmetry we get 2 unique singlet orbits

        atoms = bulk("Al", 'sc', a=1).repeat(4)
        atoms.pbc = [True, True, False]

        orbit_list = OrbitList(atoms, [0])
        cluster_space = ClusterSpace(atoms, [0], ["Al", "W"])
        self.assertEqual(len(orbit_list), 2)

        # Pbc in only one direction, so we effectively have
        # a rod-particle the unique sites are only at the surface
        # so 1 corner, one side and one bulk
        atoms = bulk("Al", 'sc', a=1).repeat(4)
        atoms.pbc = [True, False, False]

        orbit_list = OrbitList(atoms, [0])
        cluster_space = ClusterSpace(atoms, [0], ["Al", "W"])
        self.assertEqual(len(orbit_list), 3)

        # Increasing size to 5 will give
        #  one corner
        # two on the side
        # the center will be 3x3 so 3 sites from that
        atoms = bulk("Al", 'sc', a=1).repeat(5)
        atoms.pbc = [True, False, False]

        orbit_list = OrbitList(atoms, [0])
        cluster_space = ClusterSpace(atoms, [0], ["Al", "W"])
        self.assertEqual(len(orbit_list), 6)

        # Test pbc = True
        atoms = bulk("Al", 'sc', a=1).repeat(3)
        atoms.pbc = True

        orbit_list = OrbitList(atoms, [0])
        cluster_space = ClusterSpace(atoms, [0], ["Al", "W"])
        self.assertEqual(len(orbit_list), 1)


if __name__ == '__main__':
    unittest.main()
