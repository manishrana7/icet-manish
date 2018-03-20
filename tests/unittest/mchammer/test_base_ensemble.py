
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

            def __init__(self, atoms, cluster_expansion, name):
                super().__init__(atoms, cluster_expansion, name)

            def do_trial_move(self):
                pass

            def _run(self):
                pass
        self.ensemble = ConcreteEnsemble(self.atoms, self.ce, 'test-ensemble')

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

    def test_test_boltzmann_constant(self):
        """
        Test boltzmann constant.
        """
        self.assertEqual(self.ensemble.kB, 8.6173303e-5)

    def test_property_cluster_expansion(self):
        """
        Test the cluster expansion property.
        """
        self.assertEqual(self.ensemble.cluster_expansion, self.ce)

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
