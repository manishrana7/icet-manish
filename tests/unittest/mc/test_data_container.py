import unittest
import tempfile
from collections import OrderedDict
import numpy as np
import pandas as pd
from ase.build import bulk
from mchammer import DataContainer
from mchammer.observers.base_observer import BaseObserver

# Create concrete child of BaseObserver for testing


class ConcreteObserver(BaseObserver):
    def __init__(self, interval, tag='ConcreteObserver'):
        super().__init__(interval, tag=tag, return_type=int)

    def get_observable(self, atoms):
        """ Return number of Al atoms. """
        return atoms.get_chemical_symbols().count('Al')


class TestDataContainer(unittest.TestCase):
    """Container for the tests of the class functionality"""

    def __init__(self, *args, **kwargs):
        super(TestDataContainer, self).__init__(*args, **kwargs)
        self.atoms = bulk('Al').repeat(4)

    def setUp(self):
        """Set up before each test case."""
        self.dc = DataContainer(atoms=self.atoms,
                                ensemble_name='test-ensemble',
                                random_seed=44)
        test_observer = ConcreteObserver(interval=10, tag='obs1')
        self.dc.add_observable(test_observer.tag)
        self.dc.add_parameter('temperature', 375.15)

        # test another type of atoms fails
        with self.assertRaises(Exception):
            DataContainer(atoms='something',
                          ensemble_name='test-ensemble',
                          random_seed=44)

    def test_structure(self):
        """Test that reference structure is added to DataContainer."""
        self.assertEqual(self.dc.structure, self.atoms)

    def test_add_observable(self):
        """Test that observables are added to DataContainer."""
        test_observer = ConcreteObserver(interval=20, tag='obs2')
        self.dc.add_observable(test_observer.tag)
        self.assertEqual(len(self.dc.observables), 2)

        # test no duplicates
        self.dc.add_observable(test_observer.tag)
        self.assertEqual(len(self.dc.observables), 2)

    def test_add_parameter(self):
        """Test that parameters are added to DataContainer."""
        self.dc.add_parameter('sro', -0.1)
        self.assertEqual(len(self.dc.parameters), 3)

        # add a list as parameters
        index_atoms = [i for i in range(len(self.atoms))]
        self.dc.add_parameter('index_atoms', index_atoms)
        self.assertEqual(len(self.dc.parameters), 4)

    def test_append_data(self):
        """Test that data is appended to DataContainer."""
        # list of observers for testing
        observers = [ConcreteObserver(interval=10, tag='obs1'),
                     ConcreteObserver(interval=20, tag='obs2')]
        min_interval = min([obs.interval for obs in observers])
        for mctrial in range(1, 101):
            if mctrial % min_interval == 0:
                row_data = OrderedDict()
                for obs in observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.dc.append(mctrial, row_data)

        self.assertEqual(self.dc.get_number_of_entries(), 10)

    def test_property_data(self):
        """ Test data property is a Pandas DataFrame object."""
        self.assertIsInstance(self.dc.data, pd.DataFrame)

    def test_property_parameters(self):
        """Test that added parameters has OrderedDict type."""
        target = OrderedDict([('seed', 44),
                              ('temperature', 375.15)])
        retval = self.dc.parameters
        self.assertEqual(retval, target)

    def test_property_observables(self):
        """Test that added observables has list type."""
        observables = ['obs1']
        self.assertEqual(self.dc.observables, observables)

    def test_property_metadata(self):
        """Test metadata property has string type."""
        for key in self.dc.metadata:
            self.assertIsInstance(self.dc.metadata[key], str)

    def test_get_data(self):
        """
        Test the data returned as a list of list and the options provided by
        this functions works as expected.
        """
        observers = [ConcreteObserver(interval=10, tag='obs1'),
                     ConcreteObserver(interval=20, tag='obs2')]

        target = [[10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0],
                  [None, 64.0, None, 64.0, None, 64.0, None, 64.0, None, 64.0]]

        min_interval = min([obs.interval for obs in observers])
        for mctrial in range(1, 101):
            if mctrial % min_interval == 0:
                row_data = OrderedDict()
                for obs in observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.dc.append(mctrial, row_data)

        retval = self.dc.get_data()
        self.assertListEqual(target, retval)

        target = [[10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0]]

        retval = self.dc.get_data(['mctrial', 'obs1', 'obs2'],
                                  fill_missing=True)
        self.assertListEqual(target, retval)

        target = [[40, 50, 60, 70], [64.0, 64.0, 64.0, 64.0],
                  [64.0, None, 64.0, None]]

        retval = self.dc.get_data(['mctrial', 'obs1', 'obs2'],
                                  interval=(40, 70))
        self.assertListEqual(target, retval)

        with self.assertRaises(AssertionError):
            self.dc.get_data(['temperature'])

    def test_reset(self):
        """
        Test that appended data is deleted after calling this function.
        """
        for mctrial in range(100):
            self.dc.append(mctrial, dict([('temperature', 100.0)]))
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
        """Test write and read functionalities of data container."""
        # append observations to the data container
        for mctrial in range(10, 101, 10):
            self.dc.append(mctrial, dict([('temperature', 100.0)]))

        temp_file = tempfile.NamedTemporaryFile()
        self.dc.write(temp_file.name)
        dc_read = self.dc.read(temp_file.name)
        # check atoms
        self.assertEqual(self.atoms, dc_read.structure)
        # check metadata
        self.assertDictEqual(self.dc.metadata, dc_read.metadata)
        # check parameters
        self.assertDictEqual(self.dc.parameters, dc_read.parameters)
        # check observables
        self.assertEqual(self.dc.observables, dc_read.observables)
        # check runtime data
        self.assertEqual(self.dc.get_number_of_entries(),
                         dc_read.get_number_of_entries())
        self.assertListEqual(self.dc.get_data(['temperature']),
                             dc_read.get_data(['temperature']))
        # check exception raises when file does not exist
        with self.assertRaises(Exception) as context:
            dc_read = self.dc.read("not_a_file")
        msg = 'File cannot be found'
        self.assertTrue(msg in str(context.exception))


if __name__ == '__main__':
    unittest.main()
