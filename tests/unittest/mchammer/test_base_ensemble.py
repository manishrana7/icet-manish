
import unittest
from icet import ClusterSpace
from icet import ClusterExpansion
from mchammer.ensembles.base_ensemble import BaseEnsemble
from ase.build import bulk
import numpy as np


class TestEnsemble(unittest.TestCase):
    """
    Container for tests of the class functionality
    """

    def __init__(self, *args, **kwargs):
        super(TestEnsemble, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al")
        cutoffs = [5, 5, 4]
        elements = ["Al", "Ga"]
        self.cs = ClusterSpace(self.atoms, cutoffs, elements)
        parameters = np.array([1.2 for _ in range(len(self.cs))])
        self.ce = ClusterExpansion(self.cs, parameters)

    def setUp(self):
        """
        Setup before each test.
        """
        # Create a concreate child of Ensemble for testing
        class ConcreteEnsemble(BaseEnsemble):

            def __init__(self, atoms, calculator, name=None, random_seed=None):
                super().__init__(atoms, calculator, name=name, random_seed=random_seed)

            def do_trial_move(self):
                pass

            def _run(self):
                pass
        self.ensemble = ConcreteEnsemble(
            self.atoms, calculator=None, name='test-ensemble', random_seed=42)

    def test_property_name(self):
        """
        Test name property.
        """

        self.assertEqual('test-ensemble', self.ensemble.name)

    def test_property_random_seed(self):
        """
        Test random seed property.
        """

        self.assertEqual(self.ensemble.random_seed, 42)

    def test_property_accepted_trials(self):
        """
        Test property accepted trials.
        """

        self.assertEqual(self.ensemble.accepted_trials, 0)
        self.ensemble.accepted_trials += 1
        self.assertEqual(self.ensemble.accepted_trials, 1)

    def test_property_totals_trials(self):
        """
        Test property accepted trials.
        """

        self.assertEqual(self.ensemble.total_trials, 0)
        self.ensemble.total_trials += 1
        self.assertEqual(self.ensemble.total_trials, 1)

    def test_property_boltzmann_constant(self):
        """
        Test boltzmann constant.
        """
        self.assertEqual(self.ensemble.boltzmann_constant, 8.6173303e-5)

    def test_property_calculator(self):
        """
        Test the calculator property.
        """
        pass

    def test_get_next_random_number(self):
        """
        Test the get_next_random_number method.
        """
        self.assertAlmostEqual(
            self.ensemble.next_random_number(), 0.6394267984578837)

    def test_run(self):
        """
        Test the run method.
        """
        pass

    def test_internal_run(self):
        """
        Test the _run method.
        """
        pass

    def test_property_structure(self):
        """
        Test the get current structure method.
        """
        self.assertEqual(self.ensemble.structure, self.atoms)

    def test_attach_observer(self):
        """
        Test the attach method.
        """
        pass

    def test_property_data_container(self):
        """
        Test the data container property.
        """
        pass


if __name__ == '__main__':
    unittest.main()
