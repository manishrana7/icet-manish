import unittest

from ase.build import bulk
from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator


class TestCECalculator(unittest.TestCase):
    """
    Container for tests of the class functionality.

    Todo
    ----
        * add property test to calculate local contribution when that
          method has been added as intended.

    """

    def __init__(self, *args, **kwargs):
        super(TestCECalculator, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al").repeat(3)
        self.cutoffs = [0,5]
        self.subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = [1.25] * params_len

        self.ce = ClusterExpansion(self.cs, params)

    def setUp(self):
        """Setup before each test."""
        self.calculator = ClusterExpansionCalculator(
            self.atoms, self.ce, name='Test CE calc')

    def test_property_cluster_expansion(self):
        """Test the cluster expansion property."""
        self.assertIsInstance(
            self.calculator.cluster_expansion, ClusterExpansion)

    def _______test_calculate_total(self):
        """Test calculating total property."""

        self.assertEqual(self.calculator.calculate_total(
            self.atoms.numbers), 283.0)
        self.assertEqual(self.calculator.cluster_expansion.predict(
            self.calculator.atoms), 283.0)

        # set some elements
        indices = [10, 2, 4, 2]
        elements = [32] * 4
        self.calculator.update_occupations(indices, elements)
        self.assertAlmostEqual(self.calculator.calculate_total(
            self.atoms.numbers), 66.96296296)
        self.assertAlmostEqual(self.calculator.cluster_expansion.predict(
            self.calculator.atoms),  66.96296296)

    def test_calculate_local_contribution(self):
        """Test calculate local contribution."""
        indices = [i for i in range(len(self.atoms))]
        indices = [1,2,3,4]
        local_contribution = self.calculator.calculate_local_contribution(
            indices, self.atoms.numbers)
        self.assertIsInstance(local_contribution, float)

        # test local contribution by comparing with differences

        initial_value_total = self.calculator.calculate_total(
            self.atoms.numbers)
        initial_value_local = self.calculator.calculate_local_contribution(
            indices, self.atoms.numbers)
        
        current_occupations = [self.atoms.numbers[i] for i in indices]
        swapped_elements = []
        for atom in current_occupations:
            if atom == 13:
                swapped_elements.append(32)
            elif atom == 32:
                swapped_elements.append(13)
            else:
                raise Exception("Found unknown element in atoms object. {}".format(atom))
        new_occupations = self.atoms.numbers.copy()
        for index, element in zip(indices, swapped_elements):
            new_occupations[index] = element
        
        new_value_total = self.calculator.calculate_total(new_occupations)
        new_value_local = self.calculator.calculate_local_contribution(indices, new_occupations)

        total_diff = new_value_total - initial_value_total
        local_diff = new_value_local - initial_value_local
        self.assertAlmostEqual(total_diff, local_diff)

    def test_internal_calc_local_contribution(self):
        """Test the internal calc local contribution."""
        indices = [1, 2, 3]
        local_contribution = 0
        for index in indices:
            local_contribution +=\
                self.calculator._calculate_local_contribution(
                    index)
        self.assertEqual(local_contribution,
                         self.calculator.calculate_local_contribution(
                             indices, occupations=self.atoms.numbers))


if __name__ == '__main__':
    unittest.main()
