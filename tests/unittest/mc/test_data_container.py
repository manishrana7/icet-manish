#!/usr/bin/env python3
import unittest
import tempfile

from datetime import datetime
from collections import OrderedDict
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
        self.dc.add_observable('energy', float)
        self.dc.add_observable('temperature', float)
        self.assertEqual(len(self.dc.observables), 2)

    def test_add_parameter(self):
        """Test that parameters are added to DataContainer."""
        self.dc.add_parameter('temperature', 273.15)
        self.dc.add_parameter('chemical-potential-difference', -0.5)
        self.assertEqual(len(self.dc.parameters), 3)
        index_atoms = [i for i in range(len(self.atoms))]
        self.dc.add_parameter('index_atoms', index_atoms)
        self.assertEqual(len(self.dc.parameters), 4)

    def test_append_data(self):
        """Test that dc is appended to DataContainer."""
        min_interval = min([obs.interval for obs in self.observers])
        for mctrial in range(1, 101):
            if mctrial % min_interval == 0:
                row_data = OrderedDict()
                for obs in self.observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.dc.append(mctrial, row_data)

        self.assertEqual(len(self.dc), 10)

    def test_parameters(self):
        """Test that added parameters has OrderedDict type."""
        target = OrderedDict([('random-seed', 44),
                              ('temperature', 375.15),
                              ('chemical-potential-difference', -0.5)])

        self.dc.add_parameter('temperature', 375.15)
        self.dc.add_parameter('chemical-potential-difference', -0.5)

        retval = self.dc.parameters
        self.assertEqual(retval, target)

    def test_observables(self):
        """Test that added observables has OrderedDict type."""
        target = OrderedDict([('temperature', float),
                              ('latt_occupation', list)])

        self.dc.add_observable('temperature', float)
        self.dc.add_observable('latt_occupation', list)
        retval = self.dc.observables

        self.assertDictEqual(retval, target)

    def test_metadata(self):
        """Test metadata."""
        for key in self.dc.metadata:
            retval = self.dc.metadata[key]
            if key in ['date-created', 'date-last-backup']:
                self.assertIsInstance(retval, datetime)
            else:
                self.assertIsInstance(retval, str)

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
        Test this functionality removes the appended dc to DataContainer.
        """
        row_data = OrderedDict()
        row_data['temperature'] = 100.0
        for mctrial in range(10, 100, 10):
            self.dc.append(mctrial, row_data)
        self.dc.reset()
        self.assertEqual(len(self.dc), 0)

    def test_len(self):
        """Test total number of rows is returned."""
        row_data = OrderedDict()
        row_data['temperature'] = 100.0
        for mctrial in range(10, 101, 10):
            self.dc.append(mctrial, row_data)
        self.assertEqual(len(self.dc), 10)

    def test_get_average(self):
        """ Test functionality. """
        pass

    def test_restart(self):
        """ Test that BaseEnsemble object is restarted from file. """
        # add an observable
        self.dc.add_observable('temperature', float)
        # append observations to the data container
        row_data = {}
        row_data['temperature'] = 100.0
        for mctrial in range(10, 101, 10):
            self.dc.append(mctrial, row_data)
        # save to file
        temp_file = tempfile.NamedTemporaryFile()
        self.dc.write(temp_file.name)

        with self.assertRaises(NotImplementedError):
            self.dc.restart(temp_file.name)

    def test_read_write_data(self):
        """ Test write and read functionalities of data container. """
        # add an observable
        self.dc.add_observable('temperature', float)
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
        self.assertDictEqual(self.dc.observables, dc_read.observables)
        # chech runtime data
        self.assertEqual(len(self.dc), len(dc_read))
        self.assertListEqual(self.dc.get_data(['mctrial', 'temperature']),
                             dc_read.get_data(['mctrial', 'temperature']))


if __name__ == '__main__':
    unittest.main()
