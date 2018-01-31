#!/usr/bin/env Python3


import unittest

from icet.core_py.orbit_list import OrbitList
from icet.core_py.permutation_matrix import PermutationMatrix


class TestOrbitList(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        from ase.build import bulk
        super(TestOrbitList, self).__init__(*args, **kwargs)
        self.cutoffs = [4.0] * 3
        self.atoms = bulk('Ag', a=4.09)

    def setUp(self):
        '''
        Instantiate class before each test.
        '''
        self.orbit_list = OrbitList(self.atoms, self.cutoffs)

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

    def test_get_primitive_structure(self):
        '''
        Testing get_orbit_list_info functionality
        '''
        self.orbit_list.get_primitive_structure()

    def test_property_permutation_matrix(self):
        '''
        Test the permutation matrix property.
        '''
        self.assertIsInstance(
            self.orbit_list.permutation_matrix, PermutationMatrix)
        self.assertEqual(
            self.orbit_list.permutation_matrix.cutoff, max(self.cutoffs))


if __name__ == '__main__':
    unittest.main()
