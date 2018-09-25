import unittest

from ase.build import bulk
from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator
from _icet import _ClusterExpansionCalculator
from icet import Structure


class TestCECalculatorBinary(unittest.TestCase):
    """
    Container for tests of the class functionality.

    Todo
    ----
        * add property test to calculate local contribution when that
          method has been added as intended.

    """

    def __init__(self, *args, **kwargs):
        super(TestCECalculatorBinary, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al", 'hcp', a=4.0, c=3.1).repeat(2)
        # self.atoms = bulk("Al",'bcc',a=4.0).repeat(3)
        # self.atoms = bulk("Al",'diamond',a=4.0).repeat(3)
        print("atoms len ", len(self.atoms))
        self.cutoffs = [4]  # [2.9]
        self.subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = [1.0] * params_len

        self.ce = ClusterExpansion(self.cs, params)

    def setUp(self):
        """Setup before each test."""
        self.atoms = bulk("Al", 'hcp', a=4.0, c=3.1).repeat(2)
        # self.atoms = bulk("Al",'diamond',a=4.0).repeat(3)
        # self.atoms = bulk("Al",'bcc',a=4.0).repeat(3)

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

    def _test_flip_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            indices = [i]
            local_diff, total_diff = self._get_energy_diffs_local_and_total(
                indices)
            self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def _test_swap_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                if j <= i:
                    continue
                indices = [i, j]
                local_diff, total_diff = \
                    self._get_energy_diffs_local_and_total(indices)
                self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def test_local_contribution_flip(self):
        """ Test potential differences when flipping."""

        # Test original occupations
        self._test_flip_changes("original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_flip_changes("Checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_flip_changes("Segregated")

    def test_local_contribution_swap(self):
        """ test correct differences when swapping. """
        # Test original occupations
        self._test_swap_changes("Original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_swap_changes("checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_swap_changes("segregated")

    def _get_energy_diffs_local_and_total(self, indices):
        """ Get energy diffs using local and total"""

        # Original occupations
        original_occupations = self.atoms.numbers.copy()
        # Initial value total energy
        initial_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers())
        # Initial value local energy
        initial_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=self.atoms.get_atomic_numbers())
        # Flip indices

        for index in indices:
            if self.atoms[index].number == 13:
                self.atoms[index].number = 32
            elif self.atoms[index].number == 32:
                self.atoms[index].number = 13
        # Calculate new total energy
        new_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers().copy())

        # Calculate new local energy
        new_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices,
            occupations=self.atoms.get_atomic_numbers().copy())

        # difference in energy according to total energy
        total_diff = new_value_total - initial_value_total

        # Difference in energy according to local energy
        local_diff = new_value_local - initial_value_local

        # Reset occupations
        self.atoms.set_atomic_numbers(original_occupations.copy())

        return local_diff, total_diff

    def test_calculate_local_contribution(self):
        """Test calculate local contribution."""
        # indices = [i for i in range(len(self.atoms))]
        indices = [3, 5]
        local_contribution = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=self.atoms.get_atomic_numbers())
        self.assertIsInstance(local_contribution, float)

        # test local contribution by comparing with differences
        original_occupations = self.atoms.numbers.copy()
        initial_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers())

        self.atoms.set_atomic_numbers(original_occupations.copy())
        initial_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=self.atoms.get_atomic_numbers())
        self.atoms.set_atomic_numbers(original_occupations.copy())
        current_occupations = [self.atoms.get_atomic_numbers()[i]
                               for i in indices]
        self.atoms.set_atomic_numbers(original_occupations.copy())
        swapped_elements = []
        for atom in current_occupations:
            if atom == 13:
                swapped_elements.append(32)
            elif atom == 32:
                swapped_elements.append(13)
            else:
                raise Exception(
                    "Found unknown element in atoms object. {}".format(atom))

        new_occupations = self.atoms.get_atomic_numbers().copy()
        for index, element in zip(indices, swapped_elements):
            new_occupations[index] = element
        self.atoms.set_atomic_numbers(new_occupations.copy())
        new_value_total = self.calculator.calculate_total(
            occupations=new_occupations.copy())
        self.atoms.set_atomic_numbers(new_occupations.copy())
        new_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=new_occupations.copy())
        self.atoms.set_atomic_numbers(new_occupations.copy())

        total_diff = new_value_total - initial_value_total
        local_diff = new_value_local - initial_value_local
        print("total_diff: {}".format(total_diff))
        print("local_diff: {}".format(local_diff))
        print("total/local = {}".format(total_diff/local_diff))
        print("local/total = {}".format(local_diff/total_diff))
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
                             local_indices=indices,
                             occupations=self.atoms.get_atomic_numbers()))

    def test_get_local_cluster_vector(self):
        """ Tests the get local clustervector method."""

        cpp_calc = _ClusterExpansionCalculator(
            self.cs, Structure.from_atoms(self.atoms))

        index = 4
        local_cv_before = cpp_calc.get_local_cluster_vector(
            self.atoms.get_atomic_numbers(), index, [])

        self.atoms[index].symbol = 'Ge'

        local_cv_after = cpp_calc.get_local_cluster_vector(
            self.atoms.get_atomic_numbers(), index, [])

        print(local_cv_before-local_cv_after)


class TestCECalculatorBinaryHCPLongCutoffSmallSystem(unittest.TestCase):
    """
    Container for tests of the class functionality.

    Todo
    ----
        * add property test to calculate local contribution when that
          method has been added as intended.

    """

    def __init__(self, *args, **kwargs):
        super(TestCECalculatorBinaryHCPLongCutoffSmallSystem,
              self).__init__(*args, **kwargs)

        self.atoms = bulk("Al", 'hcp', a=4.0, c=3.1).repeat(2)
        print("atoms len ", len(self.atoms))
        self.cutoffs = [12]  # [2.9]
        self.subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = [1.0] * params_len
        # params[0] = 0
        # params[1] = 0

        self.ce = ClusterExpansion(self.cs, params)

    def setUp(self):
        """Setup before each test."""
        self.atoms = bulk("Al", 'hcp', a=4.0, c=3.1).repeat(2)

        self.calculator = ClusterExpansionCalculator(
            self.atoms, self.ce, name='Test CE calc')

    def _test_flip_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            indices = [i]
            local_diff, total_diff = self._get_energy_diffs_local_and_total(
                indices)
            msg += ", indices " + str(indices) + ", len of atoms " + str(len(self.atoms))
            self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def _test_swap_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                if j <= i:
                    continue
                indices = [i, j]
                local_diff, total_diff = \
                    self._get_energy_diffs_local_and_total(indices)
                msg += ", indices " + str(indices) + ", len of atoms " + str(len(self.atoms))
                self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def test_local_contribution_flip(self):
        """ Test potential differences when flipping."""

        # Test original occupations
        self._test_flip_changes("original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_flip_changes("Checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_flip_changes("Segregated")

    def test_local_contribution_swap(self):
        """ test correct differences when swapping. """
        # Test original occupations
        self._test_swap_changes("Original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_swap_changes("checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_swap_changes("segregated")

    def _get_energy_diffs_local_and_total(self, indices):
        """ Get energy diffs using local and total"""

        # Original occupations
        original_occupations = self.atoms.numbers.copy()
        # Initial value total energy
        initial_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers())
        # Initial value local energy
        initial_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=self.atoms.get_atomic_numbers())
        # Flip indices

        for index in indices:
            if self.atoms[index].number == 13:
                self.atoms[index].number = 32
            elif self.atoms[index].number == 32:
                self.atoms[index].number = 13
        # Calculate new total energy
        new_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers().copy())

        # Calculate new local energy
        new_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices,
            occupations=self.atoms.get_atomic_numbers().copy())

        # difference in energy according to total energy
        total_diff = new_value_total - initial_value_total

        # Difference in energy according to local energy
        local_diff = new_value_local - initial_value_local

        # Reset occupations
        self.atoms.set_atomic_numbers(original_occupations.copy())

        return local_diff, total_diff


class TestCECalculatorBinaryBCCLongCutoffSmallSystem(unittest.TestCase):
    """
    Container for tests of the class functionality.

    Todo
    ----
        * add property test to calculate local contribution when that
          method has been added as intended.

    """

    def __init__(self, *args, **kwargs):
        super(TestCECalculatorBinaryBCCLongCutoffSmallSystem,
              self).__init__(*args, **kwargs)

        self.atoms = bulk("Al",'bcc',a=4.0).repeat(2)        
        self.cutoffs = [17]  # [2.9]
        self.subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = [1.0] * params_len

        self.ce = ClusterExpansion(self.cs, params)

    def setUp(self):
        """Setup before each test."""
        self.atoms = bulk("Al",'bcc',a=4.0).repeat(2)
        self.calculator = ClusterExpansionCalculator(
            self.atoms, self.ce, name='Test CE calc')

    def _test_flip_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            indices = [i]
            local_diff, total_diff = self._get_energy_diffs_local_and_total(
                indices)
            self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def _test_swap_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                if j <= i:
                    continue
                indices = [i, j]
                local_diff, total_diff = \
                    self._get_energy_diffs_local_and_total(indices)
                self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def test_local_contribution_flip(self):
        """ Test potential differences when flipping."""

        # Test original occupations
        self._test_flip_changes("original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_flip_changes("Checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_flip_changes("Segregated")

    def test_local_contribution_swap(self):
        """ test correct differences when swapping. """
        # Test original occupations
        self._test_swap_changes("Original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_swap_changes("checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_swap_changes("segregated")

    def _get_energy_diffs_local_and_total(self, indices):
        """ Get energy diffs using local and total"""

        # Original occupations
        original_occupations = self.atoms.numbers.copy()
        # Initial value total energy
        initial_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers())
        # Initial value local energy
        initial_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=self.atoms.get_atomic_numbers())
        # Flip indices

        for index in indices:
            if self.atoms[index].number == 13:
                self.atoms[index].number = 32
            elif self.atoms[index].number == 32:
                self.atoms[index].number = 13
        # Calculate new total energy
        new_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers().copy())

        # Calculate new local energy
        new_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices,
            occupations=self.atoms.get_atomic_numbers().copy())

        # difference in energy according to total energy
        total_diff = new_value_total - initial_value_total

        # Difference in energy according to local energy
        local_diff = new_value_local - initial_value_local

        # Reset occupations
        self.atoms.set_atomic_numbers(original_occupations.copy())

        return local_diff, total_diff



class TestCECalculatorTernaryBCCLongCutoffSmallSystem(unittest.TestCase):
    """
    Container for tests of the class functionality.

    Todo
    ----
        * add property test to calculate local contribution when that
          method has been added as intended.

    """

    def __init__(self, *args, **kwargs):
        super(TestCECalculatorTernaryBCCLongCutoffSmallSystem,
              self).__init__(*args, **kwargs)

        self.atoms = bulk("Al",'bcc',a=4.0).repeat(2)        
        self.cutoffs = [17]
        self.subelements = ['Al', 'Ge','Si']
        self.cs = ClusterSpace(self.atoms, self.cutoffs, self.subelements)
        params_len = self.cs.get_cluster_space_size()
        params = [1.0] * params_len

        self.ce = ClusterExpansion(self.cs, params)

    def setUp(self):
        """Setup before each test."""
        self.atoms = bulk("Al",'bcc',a=4.0).repeat(2)
        self.calculator = ClusterExpansionCalculator(
            self.atoms, self.ce, name='Test CE calc')

    def _test_flip_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            indices = [i]
            local_diff, total_diff = self._get_energy_diffs_local_and_total(
                indices)
            self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def _test_swap_changes(self, msg):
        """ Test differences when flipping """
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                if j <= i:
                    continue
                indices = [i, j]
                local_diff, total_diff = \
                    self._get_energy_diffs_local_and_total(indices)
                self.assertAlmostEqual(total_diff, local_diff, msg=msg)

    def test_local_contribution_flip(self):
        """ Test potential differences when flipping."""

        # Test original occupations
        self._test_flip_changes("original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_flip_changes("Checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_flip_changes("Segregated")

    def test_local_contribution_swap(self):
        """ test correct differences when swapping. """
        # Test original occupations
        self._test_swap_changes("Original occupations")

        # Test checkerboard-ish
        for i in range(len(self.atoms)):
            if i % 2 == 0:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32

        self._test_swap_changes("checkerboard")

        # Test segregated-ish
        for i in range(len(self.atoms)):
            if i < len(self.atoms)/2:
                self.atoms[i].number = 13
            else:
                self.atoms[i].number = 32
        self._test_swap_changes("segregated")

    def _get_energy_diffs_local_and_total(self, indices):
        """ Get energy diffs using local and total"""

        # Original occupations
        original_occupations = self.atoms.numbers.copy()
        # Initial value total energy
        initial_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers())
        # Initial value local energy
        initial_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices, occupations=self.atoms.get_atomic_numbers())
        # Flip indices

        for index in indices:
            if self.atoms[index].number == 13:
                self.atoms[index].number = 32
            elif self.atoms[index].number == 32:
                self.atoms[index].number = 13
        # Calculate new total energy
        new_value_total = self.calculator.calculate_total(
            occupations=self.atoms.get_atomic_numbers().copy())

        # Calculate new local energy
        new_value_local = self.calculator.calculate_local_contribution(
            local_indices=indices,
            occupations=self.atoms.get_atomic_numbers().copy())

        # difference in energy according to total energy
        total_diff = new_value_total - initial_value_total

        # Difference in energy according to local energy
        local_diff = new_value_local - initial_value_local

        # Reset occupations
        self.atoms.set_atomic_numbers(original_occupations.copy())

        return local_diff, total_diff



if __name__ == '__main__':
    unittest.main()
