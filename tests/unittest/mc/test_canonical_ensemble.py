import unittest

import numpy as np
from ase.build import bulk

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator
from mchammer.ensembles.canonical_ensemble import CanonicalEnsemble
from test_base_ensemble import ParakeetObserver


class TestEnsemble(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestEnsemble, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al").repeat(3)
        for i, atom in enumerate(self.atoms):
            if i % 2 == 0:
                atom.symbol = "Ga"
        cutoffs = [5, 5, 4]
        elements = ["Al", "Ga"]
        self.cs = ClusterSpace(self.atoms, cutoffs, elements)
        parameters = np.array([1.2 for _ in range(len(self.cs))])
        self.ce = ClusterExpansion(self.cs, parameters)
        self.temperature = 100.0

    def setUp(self):
        """Setup before each test."""
        self.calculator = ClusterExpansionCalculator(self.atoms, self.ce)
        self.ensemble = CanonicalEnsemble(
            calculator=self.calculator, atoms=self.atoms, name='test-ensemble',
            random_seed=42, temperature=self.temperature)

        # Create an observer for testing.
        observer = ParakeetObserver(interval=7)
        self.ensemble.attach_observer(observer)
        observer = ParakeetObserver(interval=14, tag='Parakeet2')
        self.ensemble.attach_observer(observer)

    def test_temperature_attribute(self):
        """Test temperature attribute."""

        self.assertEqual(self.ensemble.temperature, self.temperature)
        self.ensemble.temperature = 300
        self.assertEqual(self.ensemble.temperature, 300)

    def test_do_trial_step(self):
        """Test the do trial step."""
        self.ensemble.do_trial_step()

        self.assertEqual(self.ensemble.total_trials, 1)

    def test_acceptance_condition(self):
        """ Test the acceptance condition method."""

        self.assertTrue(self.ensemble._acceptance_condition(-10.0))

        # at least run it for positive energy diff
        self.ensemble._acceptance_condition(10.0)

    def test_init_without_temperature(self):
        """ Test init without temperature."""
        with self.assertRaises(KeyError) as context:
            ensemble = CanonicalEnsemble(
                calculator=self.calculator, atoms=self.atoms, name='test-ensemble',
                random_seed=42)
        self.assertTrue(
            "Temperature needs to be set in canonical ensemble" in str(context.exception))


if __name__ == '__main__':
    unittest.main()
