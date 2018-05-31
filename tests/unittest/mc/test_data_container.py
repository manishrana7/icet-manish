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
        super().__init__(interval, return_type=int, tag=tag)

    def get_observable(self, atoms):
        """ Return number of Al atoms. """
        return atoms.get_chemical_symbols().count('Al')


class TestDataContainer(unittest.TestCase):
    """Container for the tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestDataContainer, self).__init__(*args, **kwargs)
        self.atoms = bulk('Al').repeat(4)

    def setUp(self):
        """Set up before each test case."""
        self.dc = DataContainer(self.atoms,
                                ensemble_name='test-ensemble',
                                random_seed=44)
        test_observer = ConcreteObserver(interval=10, tag='obs1')
        self.dc.add_observable(test_observer.tag)
        self.dc.add_parameter('temperature', 375.15)

    def test_init(self):
        """Test initializing DataContainer."""
        self.assertIsInstance(self.dc, DataContainer)

        # test fails with a non ASE Atoms type
        with self.assertRaises(Exception):
            DataContainer('atoms', 'test-ensemble', 44)

    def test_structure(self):
        """Test reference structure property."""
        self.assertEqual(self.dc.structure, self.atoms)

    def test_add_observable(self):
        """Test add observable functionality."""
        test_observer = ConcreteObserver(interval=20, tag='obs2')
        self.dc.add_observable(test_observer.tag)
        self.assertEqual(len(self.dc.observables), 2)

        # test no duplicates
        self.dc.add_observable(test_observer.tag)
        self.assertEqual(len(self.dc.observables), 2)

    def test_add_parameter(self):
        """Test add parameter functionality."""
        self.dc.add_parameter('sro', -0.1)

        # add a list as parameters
        index_atoms = [i for i in range(len(self.atoms))]
        self.dc.add_parameter('index_atoms', index_atoms)
        self.assertEqual(len(self.dc.parameters), 4)

    def test_append_data(self):
        """Test append data functionality."""
        # list of observers for testing
        observers = [ConcreteObserver(interval=10, tag='obs1'),
                     ConcreteObserver(interval=20, tag='obs2')]
        min_interval = min([obs.interval for obs in observers])

        # append data from observers
        for mctrial in range(100):
            if mctrial % min_interval == 0:
                row_data = {}
                for obs in observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.dc.append(mctrial, row_data)

        self.assertEqual(self.dc.get_number_of_entries(), 10)

        # append list type data
        row_data = {}
        row_data['occupation_vector'] = [1, 3, 7, 11]
        self.dc.append(100, row_data)
        self.assertEqual(
            self.dc.get_number_of_entries('occupation_vector'), 1)

    def test_property_data(self):
        """ Test data property."""
        self.assertIsInstance(self.dc.data, pd.DataFrame)

    def test_property_parameters(self):
        """Test parameters property."""
        self.assertEqual(self.dc.parameters,
                         OrderedDict([('seed', 44),
                                      ('temperature', 375.15)]))

    def test_property_observables(self):
        """Test observables property."""
        self.assertListEqual(self.dc.observables, ['obs1'])

    def test_property_metadata(self):
        """Test metadata property."""
        for key in self.dc.metadata:
            self.assertIsInstance(self.dc.metadata[key], str)

    def test_get_data(self):
        """
        Test the returned data is a list of list and the options provided by
        the method works as expected.
        """
        observers = [ConcreteObserver(interval=10, tag='obs1'),
                     ConcreteObserver(interval=20, tag='obs2')]

        mctrials = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        
        # append data
        min_interval = min([obs.interval for obs in observers])
        for mctrial in mctrials:
            if mctrial % min_interval == 0:
                row_data = {}
                for obs in observers:
                    if mctrial % obs.interval == 0:
                        observable = obs.get_observable(self.atoms)
                        row_data[obs.tag] = observable
                self.dc.append(mctrial, row_data)

        obs1_vals = \
            [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0] 
        obs2_vals = \ 
            [None, 64.0, None, 64.0, None, 64.0, None, 64.0, None, 64.0]

        retval = self.dc.get_data()
        self.assertEqual(retval, (mctrials, obs1_vals, obs2_vals))

        # with filling_missing = True
        target = [[10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0],
                  [64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0, 64.0]]

        retval = self.dc.get_data(fill_missing=True)
        self.assertEqual(target, retval)

        # with a given start and stop
        target = [[40, 50, 60, 70], [64.0, 64.0, 64.0, 64.0],
                  [64.0, None, 64.0, None]]

        retval = self.dc.get_data(start=40, stop=70)
        self.assertEqual(target, retval)

        # with an interval
        target = [[10, 30, 50, 70, 90]]
        retval = self.dc.get_data(['mctrial'], interval=2)
        self.assertEqual(target, retval)

        # test fails for non-stock data
        with self.assertRaises(AssertionError):
            self.dc.get_data(['temperature'])

        # append list type data
        row_data = {}
        row_data['occupation_vector'] = [1, 3, 7, 11]
        self.dc.append(200, row_data)

        # check append data
        retval = \
            self.dc.get_data(['occupation_vector'], interval=(200, 200))[0][0]
        self.assertEqual(retval, row_data['occupation_vector'])

    def test_reset(self):
        """Test appended data is cleared."""
        # add some data first
        for mctrial in range(100):
            self.dc.append(mctrial, dict([('temperature', 100.0)]))
        # clears data
        self.dc.reset()
        self.assertEqual(self.dc.get_number_of_entries(), 0)

    def test_get_number_of_entries(self):
        """Test number of entries is returned from function."""
        row_data = [100, np.nan, 1000, np.nan]
        for mctrial, data in zip([1, 2, 3, 4], row_data):
            self.dc.append(mctrial, dict([('temperature', data)]))

        self.assertEqual(self.dc.get_number_of_entries('temperature'), 2)

        # test total number of entries
        self.assertEqual(self.dc.get_number_of_entries(), 4)

    def test_get_average(self):
        """Test get average functionality."""
        n_iter, mu, sigma = 100, 1.0, 0.1
        np.random.seed(12)
        obs_val = np.random.normal(mu, sigma, n_iter).tolist()

        # append data for testing
        for mctrial in range(n_iter):
            self.dc.append(mctrial, record={'obs1': obs_val[mctrial]})

        # get average over all mctrials
        mean, std = self.dc.get_average('obs1')
        self.assertAlmostEqual(mean, 0.9855693, places=7)
        self.assertAlmostEqual(std, 0.1051220, places=7)

        # get average over slice of data
        mean, std = self.dc.get_average('obs1', start=60)
        self.assertAlmostEqual(mean, 0.9851106, places=7)
        self.assertAlmostEqual(std, 0.0993846, places=7)

        mean, std = self.dc.get_average('obs1', stop=60)
        self.assertAlmostEqual(mean, 0.9876534, places=7)
        self.assertAlmostEqual(std, 0.1095718, places=7)

        mean, std = self.dc.get_average('obs1', start=40, stop=60)
        self.assertAlmostEqual(mean, 1.0137074, places=7)
        self.assertAlmostEqual(std, 0.1152604, places=7)

    def test_read_and_write(self):
        """Test write and read functionalities of data container."""

        # append data for testing
        self.dc.add_observable('sro')
        row_data = {}
        row_data['obs1'] = 64
        row_data['occupation_vector'] = [1, 3, 7, 11]
        for mctrial in range(1, 101):
            self.dc.append(mctrial, row_data)

        temp_file = tempfile.NamedTemporaryFile()

        # check before a non-tar file
        with self.assertRaises(ValueError) as context:
            self.dc.read(temp_file.name)
        self.assertTrue('{} is not a tar file'.format(str(temp_file.name))
                        in str(context.exception))

        # save to file
        self.dc.write(temp_file.name)

        # read from file object
        dc_read = self.dc.read(temp_file)

        # check properties and metadata
        self.assertEqual(self.atoms, dc_read.structure)
        self.assertEqual(self.dc.metadata, dc_read.metadata)
        self.assertEqual(self.dc.parameters, dc_read.parameters)
        self.assertEqual(self.dc.observables, dc_read.observables)

        # check data
        pd.testing.assert_frame_equal(
            self.dc.data, dc_read.data, check_dtype=False)

        # check exception raises when file does not exist
        with self.assertRaises(FileNotFoundError):
            dc_read = self.dc.read("not_found")
        temp_file.close()


if __name__ == '__main__':
    unittest.main()
