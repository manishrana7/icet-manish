#!/usr/bin/env python3

import unittest

from mchammer import BaseEnsemble, DataContainer
from icet import ClusterSpace, ClusterExpansion
from ase.build import bulk
from datetime import datetime
from collections import OrderedDict


class TestDataContainer(unittest.TestCase):
    """
    Container for the tests of the class functionality

    Todo
    ----
    * Replace current "observers" with instance of real observers.
    """
    def __init__(self, *args, **kwargs):
        super(TestDataContainer, self).__init__(*args, **kwargs)
        subelements = ['Ag', 'Pd', 'Au']
        cutoffs = [1.4] * 3
        atoms_prim = bulk("Al")
        cluster_space = ClusterSpace(atoms_prim, cutoffs, subelements)
        params = list(range(len(cluster_space)))

        self.observers = ['termostat', 'other_obs']
        self.intervals = [100, 200]
        self.obs_tag = ['temperature', 'energy']
        self.obs_observable = [100.0, 100.0]

        self.ce = ClusterExpansion(cluster_space, params)
        # self.atoms = self.ce.primitive_structure.repeat(4)
        self.atoms = atoms_prim.repeat(4)

    def setUp(self):
        """Set up all test cases with this."""
        SomeEnsemble = BaseEnsemble(atoms=self.atoms,
                                    calculator=self.ce)
        self.data = DataContainer(SomeEnsemble)

    def test_add_structure(self):
        """Test that reference structure is added to DataContainer."""
        self.data.add_structure(self.atoms)
        self.assertEqual(self.data.structure, self.atoms)
        # test another type fails
        with self.assertRaises(Exception):
            self.data.add_structure('something')

    def test_add_observable(self):
        """Test that observables are added to DataContainer."""
        self.data.add_observable('energy', float)
        self.data.add_observable('temperature', float)
        self.assertEqual(len(self.data.observables), 2)

    def test_add_parameter(self):
        """Test that parameters are added to DataContainer"""
        self.data.add_parameter('temperature', 1000.0)
        self.data.add_parameter('chemical-potential-difference', -0.5)
        self.assertEqual(len(self.data.parameters), 2)

    def test_append_data(self):
        """Test that data is appended to DataContainer."""
        for mctrial in range(1, 1001):
            if mctrial % min(self.intervals) == 0:
                row_data = OrderedDict()
                for k, obs in enumerate(self.observers):
                    if mctrial % self.intervals[k] == 0:
                        tag = self.obs_tag[k]
                        observable = self.obs_observable[k]
                        row_data[tag] = observable
                self.data.append(mctrial, row_data)

        self.assertEqual(len(self.data), 10)

    def test_parameters(self):
        """Test that added parameters has OrderedDict type."""
        target = OrderedDict([('temperature', 300.0),
                              ('chemical-potential-difference', -0.5)])

        self.data.add_parameter('temperature', 300.0)
        self.data.add_parameter('chemical-potential-difference', -0.5)

        retval = self.data.parameters
        self.assertEqual(retval, target)

    def test_observables(self):
        """
        Test that added observables has OrderedDict type.

        Todo
        ----
        * Add seed or PRNG, interval, minimum interval
        """
        target = OrderedDict([('energy', float),
                              ('temperature', float),
                              ('sro1', list)])

        self.data.add_observable('energy', float)
        self.data.add_observable('temperature', float)
        self.data.add_observable('sro1', list)
        retval = self.data.observables

        self.assertEqual(retval, target)

    def test_metadata(self):
        """Test metadata."""
        for key in self.data.metadata:
            retval = self.data.metadata[key]
            if key == 'date-created':
                self.assertIsInstance(retval, datetime)
            else:
                self.assertIsInstance(retval, str)

    def test_get_data(self):
        """
        Test that data is returned from DataContainer has right format
        (list of lists).
        """
        target = [[100.0, 100.0], [200.0, 100.0], [300.0, 100.0],
                  [400.0, 100.0], [500.0, 100.0], [600.0, 100.0],
                  [700.0, 100.0], [800.0, 100.0], [900.0, 100.0],
                  [1000.0, 100.0]]

        for mctrial in range(1, 1001):
            if mctrial % min(self.intervals) == 0:
                row_data = OrderedDict()
                for k, obs in enumerate(self.observers):
                    if mctrial % self.intervals[k] == 0:
                        tag = self.obs_tag[k]
                        observable = self.obs_observable[k]
                        row_data[tag] = observable
                self.data.append(mctrial, row_data)

        retval = self.data.get_data(['mctrial', 'temperature'])
        self.assertListEqual(target, retval)

        with self.assertRaises(AssertionError):
            self.data.get_data(['concentration'])

    def test_reset(self):
        """
        Test this functionality removes the appended data to DataContainer.
        """
        row_data = OrderedDict()
        row_data['temperature'] = 100.0
        for mctrial in range(10, 100, 10):
            self.data.append(mctrial, row_data)
        self.data.reset()
        self.assertEqual(len(self.data), 0)

    def test_len(self):
        """Test total number of rows is returned."""
        row_data = OrderedDict()
        row_data['temperature'] = 100.0
        for mctrial in range(10, 101, 10):
            self.data.append(mctrial, row_data)
        self.assertEqual(len(self.data), 10)

    def test_get_average(self):
        """Test functionality."""
        pass

    def test_load_data(self):
        """Test read data container from file."""
        pass

    def test_save_data(self):
        """Test data container is saved in file."""
        pass


if __name__ == '__main__':
    unittest.main()
