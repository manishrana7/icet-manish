import unittest

from ase.build import bulk
from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator


class TestCECalculator(unittest.TestCase):
    """
    Container for tests of the class functionality.

    TODO:
        * add property test to calculate local contribution when that
          method has been added as intentional.

    """

    def __init__(self, *args, **kwargs):
        super(TestCECalculator, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al").repeat(3)
        self.cutoffs = [6, 6, 5]
        self.subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = list(range(params_len))

        self.ce = ClusterExpansion(self.cs, params)

    def setUp(self):
        """Setup before each test."""
        self.calculator = ClusterExpansionCalculator(
            self.atoms, self.ce, name='Test CE calc')

    def test_property_cluster_expansion(self):
        """
        Test the cluster expansion property.
        """
        self.assertIsInstance(
            self.calculator.cluster_expansion, ClusterExpansion)

    def test_calculate_total(self):
        """
        Test calculating total property.
        """

        self.assertEqual(self.calculator.calculate_total(), 283.0)
        self.assertEqual(self.calculator.cluster_expansion.predict(
            self.calculator.atoms), 283.0)

        # set some elements
        indices = [10, 2, 4, 2]
        elements = [32] * 4
        self.calculator.set_elements(indices, elements)
        self.assertAlmostEqual(self.calculator.calculate_total(), 66.96296296)
        self.assertAlmostEqual(self.calculator.cluster_expansion.predict(
            self.calculator.atoms),  66.96296296)

    def test_calculate_local_contribution(self):
        """
        Test calculate local contribution.
        """
        indices = [1, 2, 3]
        local_contribution = self.calculator.calculate_local_contribution(
            indices)
        self.assertIsInstance(local_contribution, float)

    def test_internal_calc_local_contribution(self):
        """
        Test the internal calc local contribution.
        """
        indices = [1, 2, 3]
        local_contribution = 0
        for index in indices:
            local_contribution +=\
                self.calculator._calculate_local_contribution(
                    index)
        self.assertEqual(local_contribution,
                         self.calculator.calculate_local_contribution(indices))


if __name__ == '__main__':
    unittest.main()
