#!/usr/bin/env Python3


import unittest

from icet.core_py.orbit import Orbit
from icet.core_py.lattice_site import LatticeSite
import itertools


class TestOrbit(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        super(TestOrbit, self).__init__(*args, **kwargs)

        self.lattice_sites = []
        indices = [i for i in range(8)]
        unitcell_offsets = []
        cartesian_product_lists = [[0., 1.], [0., 1.], [0., 1.]]
        for element in itertools.product(*cartesian_product_lists):
            unitcell_offsets.append(list(element))
        self.lattice_sites = [[LatticeSite(index, unitcell_offset)]
                              for index, unitcell_offset in
                              zip(indices, unitcell_offsets)]

    def setUp(self):
        '''
        Instantiate class before each test.
        '''
        self.orbit = Orbit()

    def test_init(self):
        '''
        Just testing that the setup
        (initialization) of tested class work
        '''
        # initialize from ASE Atoms
        orbit = Orbit()
        self.assertIsInstance(orbit, Orbit)

    def test_property_equivalent_sites(self):
        '''
        Testing get_orbit_list_info functionality
        '''
        self.assertEqual(self.orbit.equivalent_sites, [])
        self.orbit.equivalent_sites = self.lattice_sites
        self.assertEqual(self.orbit.equivalent_sites, self.lattice_sites)

    def test_property_representative_cluster(self):
        """
        Test getting the representative Cluster
        TODO
        ----
        * implement this more thoroughly when cluster is added
        """
        self.orbit.representative_cluster

    def test_property_representative_sites(self):
        """
        Test getting the representative sites.
        """
        with self.assertRaises(IndexError):
            self.orbit.representative_sites

        self.orbit.equivalent_sites = self.lattice_sites
        self.assertEqual(self.orbit.representative_sites,
                         self.lattice_sites[0])

    def test_property_order(self):
        """
        Test getting the order from an orbit.
        """
        with self.assertRaises(IndexError):
            self.orbit.order

        self.orbit.equivalent_sites = self.lattice_sites
        self.assertEqual(self.orbit.order,
                         len(self.lattice_sites[0]))

    def test_len(self):
        """
        Test len of orbit.
        """
        self.assertEqual(len(self.orbit), 0)
        self.orbit.equivalent_sites = self.lattice_sites
        self.assertEqual(len(self.orbit), len(self.lattice_sites))

    def test_property_geometrical_size(self):
        """
        Test the property geometrical size.
        """
        self.orbit.geometrical_size

    def test_sort(self):
        """
        Test the sort functionality.
        """
        self.orbit.equivalent_sites = sorted(self.lattice_sites,
                                             reverse=True)
        self.orbit.sort()
        self.assertEqual(self.orbit.equivalent_sites,
                         sorted(self.lattice_sites))


if __name__ == '__main__':
    unittest.main()
