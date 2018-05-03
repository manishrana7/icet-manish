import unittest

import numpy as np
from ase.build import bulk

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator
from mchammer.ensembles.base_ensemble import BaseEnsemble
from mchammer.observers.base_observer import BaseObserver


class ParakeetObserver(BaseObserver):
    """Parakeet says 2.63323e+20."""

    def __init__(self, interval, tag='Parakeet'):
        super().__init__(interval=interval, return_type=float, tag=tag)

    def get_observable(self, atoms):  # noqa
        """Say 2.63323e+20."""
        return 2.63323e+20

# Create a concrete child of Ensemble for testing


class ConcreteEnsemble(BaseEnsemble):

    def __init__(self, calculator, atoms, name=None, random_seed=None):
        super().__init__(calculator, atoms=atoms, name=name,
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

        # runt it again to check that step is the same
        n_iters = 50
        self.ensemble.run(n_iters, reset_step=True)
        self.assertEqual(self.ensemble.step, 50)

        # runt it yet again to check that step accumulates
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
            ParakeetObserver(interval=15, tag='test_Parakeet'))
        self.assertEqual(len(self.ensemble.observers), 3)
        self.assertEqual(self.ensemble.observers['test_Parakeet'].interval, 15)
        self.assertEqual(
            self.ensemble.observers['test_Parakeet'].tag, 'test_Parakeet')

    def test_property_data_container(self):
        """Test the data container property."""
        pass

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

    def test_init_without_atoms(self):
        """Test ensemble init without atoms parameter."""

        with self.assertRaises(Exception) as context:
            ensemble = ConcreteEnsemble(  # noqa
            calculator=self.calculator, atoms=None,
            name='test-ensemble', random_seed=42)

        self.assertTrue("atoms need to be set"
                        in str(context.exception))

    def test_get_property_change(self):
        """Test the get property change method."""

        initial_occupations = self.ensemble.configuration.occupations

        indices = [0, 1, 2, 3, 4]
        elements = [13, 31, 13, 31, 13]

        prop_diff = self.ensemble.get_property_change(indices, elements)
        self.assertAlmostEqual(prop_diff, 10.37037037037037)

        # Test that the method doesn't change the occupation.
        self.assertListEqual(list(initial_occupations),
                             list(self.ensemble.configuration.occupations))


if __name__ == '__main__':
    unittest.main()
