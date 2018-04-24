#!/usr/bin/env python3
import unittest
import tempfile
from collections import OrderedDict
import numpy as np
from ase.build import bulk
from mchammer import DataContainer
from mchammer.observers.base_observer import BaseObserver


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
        self.dc = DataContainer(atoms=self.atoms,
                                ensemble_name='some-ensemble',
                                random_seed=44)
        # test another type of atoms fails
        with self.assertRaises(Exception):
            self.dc = DataContainer(atoms='something',
                                    ensemble_name='some-ensemble',
                                    random_seed=44)

    def test_structure(self):
        """Test that reference structure is added to DataContainer."""
        self.assertEqual(self.dc.structure, self.atoms)

    def test_add_observable(self):
        """Test that observables are added to DataContainer."""
        self.dc.add_observable('energy')
        self.dc.add_observable('temperature')
        self.assertEqual(len(self.dc.observables), 2)

    def test_add_parameter(self):
        """Test that parameters are added to DataContainer."""
        self.dc.add_parameter('temperature', 273.15)
        self.dc.add_parameter('chemical_potential_difference', -0.5)
        self.assertEqual(len(self.dc.parameters), 3)
        index_atoms = [i for i in range(len(self.atoms))]
        self.dc.add_parameter('index_atoms', index_atoms)
        self.assertEqual(len(self.dc.parameters), 4)

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
                self.dc.append(mctrial, row_data)

        self.assertEqual(self.dc.get_number_of_entries(), 10)

    def test_parameters(self):
        """Test that added parameters has OrderedDict type."""
        target = OrderedDict([('seed', 44),
                              ('temperature', 375.15),
                              ('sro', -0.1)])

        self.dc.add_parameter('temperature', 375.15)
        self.dc.add_parameter('sro', -0.1)

        retval = self.dc.parameters
        self.assertEqual(retval, target)

    def test_observables(self):
        """Test that added observables has OrderedDict type."""
        observables = ['temperature', 'sro']
        for obs in observables:
            self.dc.add_observable(obs)
        self.assertEqual(self.dc.observables, observables)

    def test_metadata(self):
        """Test metadata."""
        for key in self.dc.metadata:
            self.assertIsInstance(self.dc.metadata[key], str)

    def test_get_data(self):
        """
        Test that dc is returned from DataContainer has right format
        (list of lists).
        """

        target = [[20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0],
                  [64.0, None, 64.0, None, 64.0, None, 64.0, None, 64.0]]

        min_interval = min([obs.interval for obs in self.observers])
        for mctrial in range(20, 101):
            if mctrial % min_interval == 0:
                row_data = OrderedDict()
                for obs in self.observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.dc.append(mctrial, row_data)

        retval = self.dc.get_data(['mctrial', 'test_observer',
                                   'other_observer'])
        print(retval)
        self.assertListEqual(target, retval)

        target = [[20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0]]

        retval = self.dc.get_data(['mctrial', 'test_observer',
                                   'other_observer'], fill_missing=True)
        self.assertListEqual(target, retval)

        target = [[40.0, 50.0, 60.0, 70.0], [64.0, 64.0, 64.0, 64.0],
                  [64.0, None, 64.0, None]]

        retval = self.dc.get_data(['mctrial', 'test_observer',
                                   'other_observer'], interval=(40, 70))
        self.assertListEqual(target, retval)

        with self.assertRaises(AssertionError):
            self.dc.get_data(['temperature'])

    def test_reset(self):
        """
        Test that appended data is deleted from DataContainer.
        """
        row_data = {}
        row_data['temperature'] = 273.15
        for mctrial in range(1, 5):
            self.dc.append(mctrial, row_data)
        self.dc.reset()
        self.assertEqual(self.dc.get_number_of_entries(), 0)

    def test_get_number_of_entries(self):
        """ Test number of entries is returned from function. """
        row_data = [100, np.nan, 1000, np.nan]
        for mctrial, data in zip(range(1, 5), row_data):
            self.dc.append(mctrial, dict([('temperature', data)]))

        self.assertEqual(self.dc.get_number_of_entries('temperature'), 2)
        self.assertEqual(self.dc.get_number_of_entries(), 4)

    def test_get_average(self):
        """ Test functionality. """
        pass

    def test_read_write_data(self):
        """ Test write and read functionalities of data container. """
        # add an observable
        self.dc.add_observable('temperature')
        # append observations to the data container
        row_data = {}
        row_data['temperature'] = 100.0
        for mctrial in range(10, 101, 10):
            self.dc.append(mctrial, row_data)

        temp_file = tempfile.NamedTemporaryFile()
        self.dc.write(temp_file.name)
        dc_read = self.dc.read(temp_file.name)
        # check metadata
        self.assertDictEqual(self.dc.metadata, dc_read.metadata)
        # check parameters
        self.assertDictEqual(self.dc.parameters, dc_read.parameters)
        # check observables
        self.assertEqual(self.dc.observables, dc_read.observables)
        # check runtime data
        self.assertEqual(self.dc.get_number_of_entries(),
                         dc_read.get_number_of_entries())
        self.assertListEqual(self.dc.get_data(['mctrial', 'temperature']),
                             dc_read.get_data(['mctrial', 'temperature']))


if __name__ == '__main__':
    unittest.main()
