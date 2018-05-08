import unittest

from ase.build import bulk
from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator

import random
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

        self.atoms = bulk("Al",'hcp',a=4.0,c=3.1).repeat(3)
        # self.atoms = bulk("Al",'bcc',a=4.0).repeat(3)
        # self.atoms = bulk("Al",'diamond',a=4.0).repeat(3)
        print("atoms len ",len(self.atoms))
        self.cutoffs = [5.1]#[2.9]
        self.subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = [1.0] * params_len
        params = [random.random() for _ in range(params_len)]

        params[0] = 0
        params[1] = 0
        print(self.cs)

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
            occupations=self.atoms.get_atomic_numbers()), 283.0)
        self.assertEqual(self.calculator.cluster_expansion.predict(
            self.calculator.atoms), 283.0)

        # set some elements
        indices = [10, 2, 4, 2]
        elements = [32] * 4
        self.calculator.update_occupations(indices, elements)
        self.assertAlmostEqual(self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers()), 66.96296296)
        self.assertAlmostEqual(self.calculator.cluster_expansion.predict(
            self.calculator.atoms),  66.96296296)
    
    def test_calculate_local_contribution(self):
        """Test calculate local contribution."""
        # indices = [i for i in range(len(self.atoms))]
        indices = [1,4,5]
        local_contribution = self.calculator.calculate_local_contribution(
            indices, self.atoms.get_atomic_numbers())
        self.assertIsInstance(local_contribution, float)

        # test local contribution by comparing with differences
        original_occupations = self.atoms.numbers.copy()
        initial_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers())
        self.atoms.set_atomic_numbers(original_occupations.copy())
        initial_value_local = self.calculator.calculate_local_contribution(
            indices, self.atoms.get_atomic_numbers())
        self.atoms.set_atomic_numbers(original_occupations.copy())
        current_occupations = [self.atoms.get_atomic_numbers()[i] for i in indices]
        self.atoms.set_atomic_numbers(original_occupations.copy())
        swapped_elements = []
        for atom in current_occupations:
            if atom == 13:
                swapped_elements.append(32)
            elif atom == 32:
                swapped_elements.append(13)
            else:
                raise Exception("Found unknown element in atoms object. {}".format(atom))
        
        new_occupations = self.atoms.get_atomic_numbers().copy()
        for index, element in zip(indices, swapped_elements):
            new_occupations[index] = element
        self.atoms.set_atomic_numbers(new_occupations.copy())
        new_value_total = self.calculator.calculate_total(occupations=new_occupations.copy())
        self.atoms.set_atomic_numbers(new_occupations.copy())
        new_value_local = self.calculator.calculate_local_contribution(indices, new_occupations.copy())
        self.atoms.set_atomic_numbers(new_occupations.copy())

        total_diff = new_value_total - initial_value_total
        local_diff = new_value_local - initial_value_local
        self.assertAlmostEqual(total_diff, local_diff)

    def __test_internal_calc_local_contribution(self):
        """Test the internal calc local contribution."""
        indices = [1, 2, 3]
        local_contribution = 0
        for index in indices:
            local_contribution +=\
                self.calculator._calculate_local_contribution(
                    index)
        self.assertEqual(local_contribution,
                         self.calculator.calculate_local_contribution(
                             indices, occupations=self.atoms.get_atomic_numbers()))


if __name__ == '__main__':
    unittest.main()
