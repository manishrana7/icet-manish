import unittest

import os
import tempfile
import numpy as np
from ase.build import bulk

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator
from mchammer.ensembles.base_ensemble import BaseEnsemble
from mchammer.observers.base_observer import BaseObserver
from mchammer import DataContainer


class ParakeetObserver(BaseObserver):
    """Parakeet says 2.63323e+20."""

    def __init__(self, interval, tag='Parakeet'):
        super().__init__(interval=interval, return_type=float, tag=tag)

    def get_observable(self, atoms):  # noqa
        """Say 2.63323e+20."""
        return 2.63323e+20


class DictObserver(BaseObserver):

    def __init__(self, interval, tag='Ayaymama'):
        super().__init__(interval=interval, return_type=dict, tag=tag)

    def get_observable(self, atoms):
        return {'value_1': 1.0, 'value_2': 2.0}

    def get_keys(self):
        return ['value_1', 'value_2']


# Create a concrete child of Ensemble for testing


class ConcreteEnsemble(BaseEnsemble):

    def __init__(self, calculator, atoms=None, name=None, data_container=None,
                 data_container_write_period=np.inf, random_seed=None):
        super().__init__(
            calculator, atoms=atoms, name=name, data_container=data_container,
            data_container_write_period=data_container_write_period,
            random_seed=random_seed)

    def do_trial_step(self):
        pass


class TestEnsemble(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestEnsemble, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al").repeat(3)
        cutoffs = [5, 5, 4]
        elements = ["Al", "Ga"]
        self.cs = ClusterSpace(self.atoms, cutoffs, elements)
        parameters = np.array([1.2 for _ in range(len(self.cs))])
        self.ce = ClusterExpansion(self.cs, parameters)

    def setUp(self):
        """Setup before each test."""
        self.calculator = ClusterExpansionCalculator(self.atoms, self.ce)
        self.ensemble = ConcreteEnsemble(
            calculator=self.calculator, atoms=self.atoms, name='test-ensemble',
            random_seed=42)

        # Create an observer for testing.
        observer = ParakeetObserver(interval=7)
        self.ensemble.attach_observer(observer)
        observer = ParakeetObserver(interval=14, tag='Parakeet2')
        self.ensemble.attach_observer(observer)

    def test_init(self):
        """Test exceptions are raised in initialisation."""
        # without atoms parameters
        with self.assertRaises(Exception) as context:
            ConcreteEnsemble(calculator=self.calculator, atoms=None,
                             name='test-ensemble', random_seed=42)

        self.assertTrue("Missing required keyword argument: atoms"
                        in str(context.exception))

        # without calculator
        with self.assertRaises(Exception) as context:
            ConcreteEnsemble(calculator=None, atoms=self.atoms,
                             name='test-ensemble', random_seed=42)

        self.assertTrue("Missing required keyword argument: calculator"
                        in str(context.exception))

    def test_property_name(self):
        """Test name property."""
        self.assertEqual('test-ensemble', self.ensemble.name)

    def test_property_random_seed(self):
        """Test random seed property."""
        self.assertEqual(self.ensemble.random_seed, 42)

    def test_property_accepted_trials(self):
        """Test property accepted trials."""
        self.assertEqual(self.ensemble.accepted_trials, 0)
        self.ensemble.accepted_trials += 1
        self.assertEqual(self.ensemble.accepted_trials, 1)

    def test_property_totals_trials(self):
        """Test property accepted trials."""
        self.assertEqual(self.ensemble.total_trials, 0)
        self.ensemble.total_trials += 1
        self.assertEqual(self.ensemble.total_trials, 1)

    def test_property_acceptance_ratio(self):
        """Test property acceptance ratio."""
        self.ensemble.total_trials = 30
        self.ensemble.accepted_trials = 15
        self.assertEqual(self.ensemble.acceptance_ratio, 0.5)

    def test_property_calculator(self):
        """Test the calculator property."""
        pass

    def test_get_next_random_number(self):
        """Test the get_next_random_number method."""
        self.assertAlmostEqual(
            self.ensemble.next_random_number(), 0.6394267984578837)

    def test_run(self):
        """Test the run method."""

        n_iters = 364
        self.ensemble.run(n_iters)
        self.assertEqual(self.ensemble.step, n_iters)
        dc_data = self.ensemble.data_container.get_data(tags=['Parakeet2'])

        number_of_observations = len([x for x in dc_data[0] if x is not None])
        # plus one since we also count step 0
        self.assertEqual(
            number_of_observations,
            n_iters // self.ensemble.observers['Parakeet2'].interval + 1)

        # run it again to check that step is the same
        n_iters = 50
        self.ensemble.run(n_iters, reset_step=True)
        self.assertEqual(self.ensemble.step, 50)

        # run it yet again to check that step accumulates
        n_iters = 10
        self.ensemble.run(n_iters, reset_step=False)
        self.ensemble.run(n_iters, reset_step=False)
        self.assertEqual(self.ensemble.step, 70)

        # Do a number of steps of continuous runs and see that
        # we get the expected number of parakeet observations.
        for i in range(30):
            self.ensemble.reset_data_container()
            run_iters = [1, 50, 100, 200, i]
            for n_iter in run_iters:
                self.ensemble.run(n_iter)
            total_iters = sum(run_iters)
            # Check that the number of iters are correct
            self.assertEqual(self.ensemble.step, total_iters)
            dc_data = self.ensemble.data_container.get_data(tags=['Parakeet2'])
            number_of_observations = len(
                [x for x in dc_data[0] if x is not None])
            # plus one since we also count step 0
            self.assertEqual(
                number_of_observations,
                total_iters //
                self.ensemble.observers['Parakeet2'].interval + 1)

    def test_run_with_dict_observer(self):
        """Test the run method with a dict observer."""
        observer = DictObserver(interval=28)
        self.ensemble.attach_observer(observer)

        n_iters = 364
        self.ensemble.run(n_iters)
        self.assertEqual(self.ensemble.step, n_iters)
        dc_data = \
            self.ensemble.data_container.get_data(tags=['value_1', 'value_2'])

        self.assertEqual(len(dc_data[0]), len(dc_data[1]))

        number_of_observations = len([x for x in dc_data[0] if x is not None])
        # plus one since we also count step 0
        self.assertEqual(
            number_of_observations,
            n_iters // self.ensemble.observers['Ayaymama'].interval + 1)

    def test_backup_file(self):
        """Test data is being saved and can be read by the ensemble."""
        # set-up ensemble with a non-inf write period
        ensemble = ConcreteEnsemble(calculator=self.calculator,
                                    atoms=self.atoms,
                                    data_container='mycontainer.dc',
                                    data_container_write_period=1e-4)

        # attach observer
        observer = ParakeetObserver(interval=14, tag='Parakeet2')
        ensemble.attach_observer(observer)

        # back-up data while run ensemble and then read the file
        try:
            n_iters = 364
            ensemble.run(n_iters)
            dc_read = DataContainer.read('mycontainer.dc')

        finally:
            os.remove('mycontainer.dc')

        # check data container
        dc_data = dc_read.get_data(tags=['Parakeet2'])
        self.assertEqual(
            len(dc_data[0]),
            n_iters // observer.interval + 1)

        # write data container to tempfile
        temp_container_file = tempfile.NamedTemporaryFile()
        dc_read.write(temp_container_file.name)

        # initialise a new ensemble with dc file
        ensemble_reloaded = \
            ConcreteEnsemble(calculator=self.calculator,
                             atoms=self.atoms,
                             data_container=temp_container_file.name)

        # check loaded data container of new ensemble
        data_dc_reloaded = \
            ensemble_reloaded.data_container.get_data(tags=['Parakeet2'])
        data_dc = \
            ensemble.data_container.get_data(tags=['Parakeet2'])
        self.assertEqual(data_dc[0], data_dc_reloaded[0])

    def test_internal_run(self):
        """Test the _run method."""
        pass

    def test_property_structure(self):
        """Test the get current structure method."""
        # need calculator for structure
        pass
        # self.assertEqual(self.ensemble.structure, self.atoms)

    def test_attach_observer(self):
        """Test the attach method."""
        self.assertEqual(len(self.ensemble.observers), 2)

        self.ensemble.attach_observer(
            ParakeetObserver(interval=10, tag='test_Parakeet'))
        self.assertEqual(self.ensemble.observers['test_Parakeet'].interval, 10)
        self.assertEqual(
            self.ensemble.observers['test_Parakeet'].tag, 'test_Parakeet')
        self.assertEqual(len(self.ensemble.observers), 3)

        # test no duplicates, this should overwrite the last Parakeet
        self.ensemble.attach_observer(
            ParakeetObserver(interval=15), tag='test_Parakeet')
        self.assertEqual(len(self.ensemble.observers), 3)
        self.assertEqual(self.ensemble.observers['test_Parakeet'].interval, 15)
        self.assertEqual(
            self.ensemble.observers['test_Parakeet'].tag, 'test_Parakeet')

    def test_property_data_container(self):
        """Test the data container property."""
        self.assertIsInstance(self.ensemble.data_container, DataContainer)

    def test_find_minimum_observation_interval(self):
        """Test the method to find the minimum observation interval."""
        pass

    def test_property_minimum_observation_interval(self):
        """Test property minimum observation interval."""
        pass

    def test_get_gcd(self):
        """Test the get gcd method."""
        input = [2, 4, 6, 8]
        target = 2
        self.assertEqual(self.ensemble._get_gcd(input), target)

        input = [20, 15, 10]
        target = 5
        self.assertEqual(self.ensemble._get_gcd(input), target)

    def test_get_property_change(self):
        """Test the get property change method."""

        initial_occupations = self.ensemble.configuration.occupations

        indices = [0, 1, 2, 3, 4]
        elements = [13, 31, 13, 31, 13]

        prop_diff = self.ensemble.get_property_change(indices, elements)
        self.assertAlmostEqual(prop_diff, 2.0740740740740735)

        # Test that the method doesn't change the occupation.
        self.assertListEqual(list(initial_occupations),
                             list(self.ensemble.configuration.occupations))

        with self.assertRaises(ValueError) as context:
            self.ensemble.update_occupations(indices, elements+[31])

        self.assertTrue(
            "List of sites and list of elements are not the same size."
            in str(context.exception))


if __name__ == '__main__':
    unittest.main()
