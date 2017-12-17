#!/usr/bin/env Python3

'''
This file contains unit tests and other tests. It can be executed by
simply executing this file from a shell prompt:

    $ ./test_cluster_space.py

In which case it will use the system's default Python version. If a specific
Python version should be used, run that Python version with this file as input,
e.g.:

    python3 test_cluster_space.py

For a description of the Python unit testing framework, see this link:
https://docs.python.org/3/library/unittest.html

When executing this file doc testing is also performed on all doc tests in
the cluster_space.py file

'''

import unittest

from icetdev import ClusterSpace
from icetdev.cluster_space import (get_singlet_info,
                                   get_singlet_configuration,
                                   get_Mi_from_dict)
from icetdev.lattice_site import LatticeSite
from ase.build import bulk
from ase import Atoms
from collections import OrderedDict


def strip_surrounding_spaces(inp):
    '''
    Helper function that removes both leading and trailing spaces from a
    multi-line string.
    '''
    from io import StringIO
    s = []
    for line in StringIO(inp):
        if len(line.strip()) == 0:
            continue
        s += [line.strip()]
    return '\n'.join(s)


def _assertEqualComplexList(self, retval, target):
    '''
    Helper function that conducts a systematic comparison of a nested list
    with dictionaries.
    '''
    self.assertIsInstance(retval, type(target))
    for row_retval, row_target in zip(retval, target):
        self.assertIsInstance(row_retval, type(row_target))
        for key, val in row_target.items():
            self.assertIn(key, row_retval)
            s = ['key: {}'.format(key)]
            s += ['type: {}'.format(type(key))]
            s += ['retval: {}'.format(row_retval[key])]
            s += ['target: {}'.format(val)]
            info = '   '.join(s)
            if isinstance(val, float):
                self.assertAlmostEqual(val, row_retval[key], places=9,
                                       msg=info)
            else:
                self.assertEqual(row_retval[key], val, msg=info)

unittest.TestCase.assertEqualComplexList = _assertEqualComplexList


class TestClusterSpace(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def setUp(self):
        '''
        Instantiate class before each test.

        '''
        self.cs = ClusterSpace(atoms_prim, cutoffs, subelements)

    def test_init(self):
        '''
        Just testing that the setup
        (initialization) of tested class work

        '''
        self.cs = ClusterSpace(atoms_prim, cutoffs, subelements)

    def test_len(self):
        '''
        Testing len functionality

        '''
        number_orbits = self.cs.__len__()
        self.assertEqual(number_orbits, len(self.cs.get_orbit_list()))

    def test_get_data(self):
        '''
        Testing get_data functionality
        '''
        parameters = list(range(len(self.cs)+1))
        retval = self.cs.get_data(parameters=parameters)
        target = [OrderedDict([('index', 0), ('order', 1),
                               ('size', 0.0),
                               ('multiplicity', 1),
                               ('orbit_index', 0),
                               ('mc_vector', [0]),
                               ('ECI', 1)]),
                  OrderedDict([('index', 1), ('order', 2),
                               ('size', 1.4460333675264896),
                               ('multiplicity', 6),
                               ('orbit_index', 1),
                               ('mc_vector', [0, 0]),
                               ('ECI', 2)]),
                  OrderedDict([('index', 2), ('order', 3),
                               ('size', 1.6697355079971996),
                               ('multiplicity', 8),
                               ('orbit_index', 2),
                               ('mc_vector', [0, 0, 0]),
                               ('ECI', 3)]),
                  OrderedDict([('index', 3), ('order', 4),
                               ('size', 1.771021950739177),
                               ('multiplicity', 2),
                               ('orbit_index', 3),
                               ('mc_vector', [0, 0, 0, 0]),
                               ('ECI', 4)])]
        self.assertEqualComplexList(retval, target)

    def test_repr(self):
        '''
        Testing repr functionality
        '''
        retval = self.cs.__repr__()
        target = '''
------------------------- Cluster Space -------------------------
 subelements: Ag Au
 cutoffs: 4.0 4.0 4.0
 number of orbits: 4
-----------------------------------------------------------------
order |  radius  | multiplicity | index | orbit |    MC vector
-----------------------------------------------------------------
  1   |   0.0000 |        1     |    0  |    0  |    [0]
  2   |   1.4460 |        6     |    1  |    1  |  [0, 0]
  3   |   1.6697 |        8     |    2  |    2  | [0, 0, 0]
  4   |   1.7710 |        2     |    3  |    3  | [0, 0, 0, 0]
-----------------------------------------------------------------
'''
        self.assertEqual(strip_surrounding_spaces(target),
                         strip_surrounding_spaces(retval))

    def test_get_cluster_vector(self):
        '''
        Testing get_cluster_vector functionality
        '''
        for atoms, target in zip(list_atoms, target_cluster_vectors):
            retval = self.cs.get_cluster_vector(atoms)
            self.assertAlmostEqual(retval, target, places=9)

    def test_get_singlet_info(self):
        '''
        Testing get_singlet_info functionality
        '''
        retval = get_singlet_info(atoms)
        target = [{'orbit_index': 0,
                   'sites': [[LatticeSite(0, [0., 0., 0.])]],
                   'multiplicity': 1,
                   'representative_site': [LatticeSite(0, [0., 0., 0.])]}]
        self.assertEqualComplexList(retval, target)
        retval1, retval2 = get_singlet_info(atoms, return_cluster_space=True)
        self.assertEqualComplexList(retval1, target)
        self.assertIsInstance(retval2, type(self.cs))

    def test_get_singlet_configuration(self):
        '''
        Testing get_singlet_configuration functionality
        '''
        retval = get_singlet_configuration(atoms_prim)
        self.assertIsInstance(retval, Atoms)
        retval = get_singlet_configuration(list_atoms[0], to_primitive=True)
        self.assertIsInstance(retval, Atoms)
        self.assertEqual(len(retval), len(atoms_prim))

    def test_get_Mi_from_dict(self):
        '''
        Testing get_Mi_from_dict functionality
        Todo
        ----
        implement a test for get_Mi_from_dict
        '''


def suite():
    test_classes_to_run = [TestClusterSpace]
    suites_list = []
    for test_class in test_classes_to_run:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(test_class)
        suites_list.append(suite)
    test_suite = unittest.TestSuite(suites_list)
    return test_suite


subelements = ['Ag', 'Au']
cutoffs = [4.0] * 3
atoms_prim = bulk('Ag')
list_atoms = []
for k in range(4):
    atoms = atoms_prim.repeat(2)
    symbols = ['Ag'] * len(atoms)
    symbols[:k] = ['Au'] * k
    atoms.set_chemical_symbols(symbols)
    list_atoms.append(atoms)
target_cluster_vectors = [
    [1.0, -1.0, 1.0, -1.0, 1.0],
    [1.0, -0.75, 0.5, -0.25, 0.0],
    [1.0, -0.5, 0.16666666666666666, 0.0, 0.0],
    [1.0, -0.25, 0.0, 0.0, 0.0]]


if __name__ == '__main__':
    unittest.main()
