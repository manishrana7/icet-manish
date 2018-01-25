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
        from ase.build import bulk
        super(TestOrbit, self).__init__(*args, **kwargs)

        self.lattice_sites = []
        indices = [i for i in range(8)]
        unitcell_offsets = []
        cartesian_product_lists = [[0., 1.], [0., 1.], [0., 1.]]
        for element in itertools.product(*cartesian_product_lists):
            unitcell_offsets.append(list(element))
        self.lattice_sites = [LatticeSite(index, unitcell_offset)
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


if __name__ == '__main__':
    unittest.main()
