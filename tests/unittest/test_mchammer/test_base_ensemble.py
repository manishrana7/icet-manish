import unittest

import os
import tempfile
import numpy as np
from ase.build import bulk
from pandas.testing import assert_frame_equal

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator
from mchammer.ensembles.base_ensemble import BaseEnsemble
from mchammer.observers.base_observer import BaseObserver
from mchammer import DataContainer


class AppleObserver(BaseObserver):
    """Apple says 2.63323e+20."""

    def __init__(self, interval, tag='Apple'):
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

    def __init__(self, atoms, calculator, name=None, data_container=None,
                 data_container_write_period=np.inf, random_seed=None,
                 ensemble_data_write_interval=None,
                 trajectory_write_interval=None):
        super().__init__(
            atoms=atoms, calculator=calculator, name=name,
            data_container=data_container,
            data_container_write_period=data_container_write_period,
            random_seed=random_seed,
            ensemble_data_write_interval=ensemble_data_write_interval,
            trajectory_write_interval=trajectory_write_interval)

    def _do_trial_step(self):
        pass


class TestEnsemble(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestEnsemble, self).__init__(*args, **kwargs)

        self.atoms = bulk('Al').repeat(3)
        cutoffs = [5, 5, 4]
        elements = ['Al', 'Ga']
        self.cs = ClusterSpace(self.atoms, cutoffs, elements)
        parameters = np.array([1.2 for _ in range(len(self.cs))])
        self.ce = ClusterExpansion(self.cs, parameters)

    def shortDescription(self):
        """Silences unittest from printing the docstrings in test cases."""
        return None

    def setUp(self):
        """Setup before each test."""
        self.calculator = ClusterExpansionCalculator(self.atoms, self.ce)
        self.ensemble = ConcreteEnsemble(
            atoms=self.atoms, calculator=self.calculator, name='test-ensemble',
            random_seed=42)

        # Create an observer for testing.
        observer = AppleObserver(interval=7)
        self.ensemble.attach_observer(observer)
        observer = AppleObserver(interval=14, tag='Apple2')
        self.ensemble.attach_observer(observer)

    def test_init(self):
        """Tests exceptions are raised in initialization."""

        with self.assertRaises(TypeError) as context:
            ConcreteEnsemble(calculator=self.calculator)
        self.assertIn("required positional argument: 'atoms'",
                      str(context.exception))

        with self.assertRaises(TypeError) as context:
            ConcreteEnsemble(atoms=self.atoms)
        self.assertIn("required positional argument: 'calculator'",
                      str(context.exception))

        # wrong path to data container file
        with self.assertRaises(FileNotFoundError) as context:
            ConcreteEnsemble(atoms=self.atoms,
                             calculator=self.calculator,
                             data_container='path/to/nowhere/mydc')

        self.assertTrue('Path to data container file does not exist:'
                        ' path/to/nowhere' in str(context.exception))

    def test_property_name(self):
        """Tests name property."""
        self.assertEqual('test-ensemble', self.ensemble.name)

    def test_property_atoms(self):
        """Tests atoms property."""
        self.assertEqual(self.atoms, self.ensemble.atoms)

    def test_property_random_seed(self):
        """Tests random seed property."""
        self.assertEqual(self.ensemble.random_seed, 42)

    def test_property_accepted_trials(self):
        """Tests property accepted trials."""
        self.assertEqual(self.ensemble._accepted_trials, 0)
        self.ensemble._accepted_trials += 1
        self.assertEqual(self.ensemble._accepted_trials, 1)

    def test_property_totals_trials(self):
        """Tests property accepted trials."""
        self.assertEqual(self.ensemble._total_trials, 0)
        self.ensemble._total_trials += 1
        self.assertEqual(self.ensemble._total_trials, 1)

    def test_property_step(self):
        """Tests property accepted trials."""
        self.assertEqual(self.ensemble._step, 0)
        self.ensemble._step += 1
        self.assertEqual(self.ensemble._step, 1)

    def test_property_acceptance_ratio(self):
        """Tests property acceptance ratio."""
        self.ensemble._total_trials = 30
        self.ensemble._accepted_trials = 15
        self.assertEqual(self.ensemble.acceptance_ratio, 0.5)

    def test_property_calculator(self):
        """Tests the calculator property."""
        pass

    def test_get_next_random_number(self):
        """Tests the get_next_random_number method."""
        self.assertAlmostEqual(
            self.ensemble._next_random_number(), 0.6394267984578837)

    def test_run(self):
        """Tests the run method."""

        n_iters = 364
        self.ensemble.run(n_iters)
        self.assertEqual(self.ensemble._step, n_iters)
        dc_data = self.ensemble.data_container.get_data('Apple2')

        number_of_observations = len([x for x in dc_data if x is not None])
        # plus one since we also count step 0
        self.assertEqual(
            number_of_observations,
            n_iters // self.ensemble.observers['Apple2'].interval + 1)

        # run it again to check that step is the same
        n_iters = 50
        self.ensemble.run(n_iters, reset_step=True)
        self.assertEqual(self.ensemble._step, 50)

        # run it yet again to check that step accumulates
        n_iters = 10
        self.ensemble.run(n_iters, reset_step=False)
        self.ensemble.run(n_iters, reset_step=False)
        self.assertEqual(self.ensemble._step, 70)

        # Do a number of steps of continuous runs and see that
        # we get the expected number of parakeet observations.
        for i in range(30):
            self.ensemble.reset_data_container()
            run_iters = [1, 50, 100, 200, i]
            for n_iter in run_iters:
                self.ensemble.run(n_iter)
            total_iters = sum(run_iters)
            # Check that the number of iters are correct
            self.assertEqual(self.ensemble._step, total_iters)
            dc_data = self.ensemble.data_container.get_data('Apple2')
            number_of_observations = len(
                [x for x in dc_data if x is not None])
            # plus one since we also count step 0
            self.assertEqual(
                number_of_observations,
                total_iters // self.ensemble.observers['Apple2'].interval + 1)

    def test_run_with_dict_observer(self):
        """Tests the run method with a dict observer."""
        observer = DictObserver(interval=28)
        self.ensemble.attach_observer(observer)

        n_iters = 364
        self.ensemble.run(n_iters)
        self.assertEqual(self.ensemble._step, n_iters)
        dc_data = \
            self.ensemble.data_container.get_data('value_1', 'value_2')

        self.assertEqual(len(dc_data[0]), len(dc_data[1]))

        number_of_observations = len([x for x in dc_data[0] if x is not None])
        # plus one since we also count step 0
        self.assertEqual(
            number_of_observations,
            n_iters // self.ensemble.observers['Ayaymama'].interval + 1)

    def test_backup_file(self):
        """Tests data is being saved and can be read by the ensemble."""
        # set-up ensemble with a non-inf write period
        ensemble = ConcreteEnsemble(atoms=self.atoms,
                                    calculator=self.calculator,
                                    name='this-ensemble',
                                    data_container='my-datacontainer.dc',
                                    data_container_write_period=1e-2,
                                    ensemble_data_write_interval=14,
                                    trajectory_write_interval=56)

        # attach observer
        observer = AppleObserver(interval=14, tag='Apple2')
        ensemble.attach_observer(observer)

        # back-up data while run ensemble and then read the file
        try:
            n_iters = 182
            ensemble.run(n_iters)
            dc_read = DataContainer.read('my-datacontainer.dc')
        finally:
            os.remove('my-datacontainer.dc')

        # check data container
        dc_data = dc_read.get_data('Apple2')
        self.assertEqual(
            len(dc_data),
            n_iters // observer.interval + 1)

        # write data container to tempfile
        temp_container_file = tempfile.NamedTemporaryFile()
        dc_read.write(temp_container_file.name)

        # initialise a new ensemble with dc file
        ensemble_reloaded = \
            ConcreteEnsemble(atoms=self.atoms,
                             calculator=self.calculator,
                             data_container=temp_container_file.name,
                             ensemble_data_write_interval=14,
                             trajectory_write_interval=56)

        assert_frame_equal(ensemble.data_container.data,
                           ensemble_reloaded.data_container.data,
                           check_dtype=False)

        # run old and new ensemble and check both data containers are equal
        try:
            n_iters = 50
            ensemble.run(n_iters)
        finally:
            os.remove('my-datacontainer.dc')

        ensemble_reloaded.attach_observer(observer)
        ensemble_reloaded.run(n_iters)

        assert_frame_equal(ensemble.data_container.data,
                           ensemble_reloaded.data_container.data,
                           check_dtype=False)

        self.assertEqual(
            ensemble_reloaded.data_container.last_state['last_step'],
            182 + 50)

    def test_internal_run(self):
        """Tests the _run method."""
        pass

    def test_attach_observer(self):
        """Tests the attach method."""
        self.assertEqual(len(self.ensemble.observers), 2)

        self.ensemble.attach_observer(
            AppleObserver(interval=10, tag='test_Apple'))
        self.assertEqual(self.ensemble.observers['test_Apple'].interval, 10)
        self.assertEqual(
            self.ensemble.observers['test_Apple'].tag, 'test_Apple')
        self.assertEqual(len(self.ensemble.observers), 3)

        # test no duplicates, this should overwrite the last Apple
        self.ensemble.attach_observer(
            AppleObserver(interval=15), tag='test_Apple')
        self.assertEqual(len(self.ensemble.observers), 3)
        self.assertEqual(self.ensemble.observers['test_Apple'].interval, 15)
        self.assertEqual(
            self.ensemble.observers['test_Apple'].tag, 'test_Apple')

        # check that correct exceptions are raised
        with self.assertRaises(TypeError) as context:
            self.ensemble.attach_observer('xyz')
        self.assertTrue('observer has the wrong type'
                        in str(context.exception))

    def test_property_data_container(self):
        """Tests the data container property."""
        self.assertIsInstance(self.ensemble.data_container, DataContainer)

    def test_find_minimum_observation_interval(self):
        """Tests the method to find the minimum observation interval."""
        pass

    def test_property_minimum_observation_interval(self):
        """Tests property minimum observation interval."""
        pass

    def test_get_gcd(self):
        """Tests the get gcd method."""
        input = [2, 4, 6, 8]
        target = 2
        self.assertEqual(self.ensemble._get_gcd(input), target)

        input = [20, 15, 10]
        target = 5
        self.assertEqual(self.ensemble._get_gcd(input), target)

    def test_get_property_change(self):
        """Tests the get property change method."""

        initial_occupations = self.ensemble.configuration.occupations

        indices = [0, 1, 2, 3, 4]
        elements = [13, 31, 13, 31, 13]

        prop_diff = self.ensemble._get_property_change(indices, elements)
        self.assertAlmostEqual(prop_diff, 56)

        # Tests that the method doesn't change the occupation.
        self.assertListEqual(list(initial_occupations),
                             list(self.ensemble.configuration.occupations))

        with self.assertRaises(ValueError) as context:
            self.ensemble.update_occupations(indices, elements + [31])

        self.assertTrue('sites and species must have the same length.'
                        in str(context.exception))

    def test__get_ensemble_data(self):
        """Tests the get ensemble data method."""
        data = self.ensemble._get_ensemble_data()

        self.assertIn('potential', data.keys())

    def test_data_container_write_period(self):
        """Tests property data container write container."""
        self.assertTrue(np.isinf(self.ensemble.data_container_write_period))

        self.ensemble.data_container_write_period = 1e-2
        self.assertAlmostEqual(self.ensemble.data_container_write_period, 1e-2)


if __name__ == '__main__':
    unittest.main()