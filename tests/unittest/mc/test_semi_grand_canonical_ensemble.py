import unittest

import numpy as np
from ase.build import bulk

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators import ClusterExpansionCalculator

from mchammer.ensembles import SemiGrandCanonicalEnsemble


class TestEnsemble(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestEnsemble, self).__init__(*args, **kwargs)

        self.atoms = bulk('Al').repeat(3)
        for i, atom in enumerate(self.atoms):
            if i % 2 == 0:
                atom.symbol = 'Ga'
        cutoffs = [5, 5, 4]
        elements = ['Al', 'Ga']
        self.chemical_potentials = {'Al': 5, 'Ga':0}
        self.cs = ClusterSpace(self.atoms, cutoffs, elements)
        parameters = parameters = np.array([1.2]*len(self.cs))
        self.ce = ClusterExpansion(self.cs, parameters)
        self.temperature = 100.0

    def setUp(self):
        """Setup before each test."""
        self.calculator = ClusterExpansionCalculator(self.atoms, self.ce)

        self.ensemble = SemiGrandCanonicalEnsemble(
            calculator=self.calculator, atoms=self.atoms, name='test-ensemble',
            random_seed=42, temperature=self.temperature, chemical_potentials=self.chemical_potentials)

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
            SemiGrandCanonicalEnsemble(
                calculator=self.calculator, atoms=self.atoms,
                name='test-ensemble', random_seed=42)
        self.assertTrue(
            "Temperature needs to be set in canonical "
            "ensemble" in str(context.exception))


if __name__ == '__main__':
    unittest.main()
