#!/usr/bin/env python3
import unittest

from mchammer import DataContainer
from mchammer.observers.base_observer import BaseObserver
from ase.build import bulk
from datetime import datetime
from collections import OrderedDict


class TestDataContainer(unittest.TestCase):
    """Container for the tests of the class functionality"""

    def __init__(self, *args, **kwargs):
        super(TestDataContainer, self).__init__(*args, **kwargs)
        self.atoms = bulk('Al').repeat(4)

        class ConcreteObserver(BaseObserver):
            def __init__(self, interval, tag='ConcreteObserver'):
                super().__init__(interval, tag=tag)

            @property
            def return_type(self):
                return int

            def get_observable(self, atoms):
                """
                Return number of Al atoms.
                """
                return atoms.get_chemical_symbols().count('Al')

        obs1 = ConcreteObserver(interval=10, tag='test_observer')
        obs2 = ConcreteObserver(interval=20, tag='other_observer')
        self.observers = [obs1, obs2]

    def setUp(self):
        """Set up before each test case."""
        self.data = DataContainer(atoms=self.atoms,
                                  ensemble_name='some-ensemble',
                                  random_seed=44)

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
        """Test that parameters are added to DataContainer."""
        self.data.add_parameter('temperature', 273.15)
        self.data.add_parameter('chemical-potential-difference', -0.5)
        self.assertEqual(len(self.data.parameters), 3)
        index_atoms = [i for i in range(len(self.atoms))]
        self.data.add_parameter('index_atoms', index_atoms)
        self.assertEqual(len(self.data.parameters), 4)

    def test_append_data(self):
        """Test that data is appended to DataContainer."""
        min_interval = min([obs.interval for obs in self.observers])
        for mctrial in range(1, 101):
            if mctrial % min_interval == 0:
                row_data = OrderedDict()
                for obs in self.observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.data.append(mctrial, row_data)

        self.assertEqual(len(self.data), 10)

    def test_parameters(self):
        """Test that added parameters has OrderedDict type."""
        target = OrderedDict([('random-seed', 44),
                              ('temperature', 375.15),
                              ('chemical-potential-difference', -0.5)])

        self.data.add_parameter('temperature', 375.15)
        self.data.add_parameter('chemical-potential-difference', -0.5)

        retval = self.data.parameters
        self.assertEqual(retval, target)

    def test_observables(self):
        """Test that added observables has OrderedDict type."""
        target = OrderedDict([('temperature', float),
                              ('latt_occupation', list)])

        self.data.add_observable('temperature', float)
        self.data.add_observable('latt_occupation', list)
        retval = self.data.observables

        self.assertDictEqual(retval, target)

    def test_get_data(self):
        """
        Test that data is returned from DataContainer has right format
        (list of lists).
        """
        target = [[20.0, 64.0, 64.0], [30.0, 64.0, None],
                  [40.0, 64.0, 64.0], [50.0, 64.0, None], [60.0, 64.0, 64.0],
                  [70.0, 64.0, None], [80.0, 64.0, 64.0], [90.0, 64.0, None],
                  [100.0, 64.0, 64.0]]

        min_interval = min([obs.interval for obs in self.observers])
        for mctrial in range(20, 101):
            if mctrial % min_interval == 0:
                row_data = OrderedDict()
                for obs in self.observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.data.append(mctrial, row_data)

        retval = self.data.get_data(['mctrial', 'test_observer',
                                     'other_observer'])
        self.assertListEqual(target, retval)

        target = [[20.0, 64.0, 64.0], [30.0, 64.0, 64.0],
                  [40.0, 64.0, 64.0], [50.0, 64.0, 64.0], [60.0, 64.0, 64.0],
                  [70.0, 64.0, 64.0], [80.0, 64.0, 64.0], [90.0, 64.0, 64.0],
                  [100.0, 64.0, 64.0]]

        retval = self.data.get_data(['mctrial', 'test_observer',
                                     'other_observer'], fill_missing=True)
        self.assertListEqual(target, retval)

        with self.assertRaises(AssertionError):
            self.data.get_data(['temperature'])

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
