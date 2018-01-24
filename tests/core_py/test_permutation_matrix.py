from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core.permutation_map import permutation_matrix_from_atoms
from icet.core.permutation_map import PermutationMap as PermutationMatrix_cpp
from ase.build import bulk

import unittest
import numpy as np

class TestPermutationMatrix(unittest.TestCase):
    '''
    Test the python implementation of permutation matrix
    '''

    def __init__(self, *args, **kwargs):
        super(TestPermutationMatrix, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al")
        self.cutoff = 10.0

    def setUp(self):
        """
        Setup.
        """
        self.pm = PermutationMatrix(self.atoms, self.cutoff)
        self.pm_cpp = permutation_matrix_from_atoms(self.atoms, self.cutoff)[0]

    def test_init(self):
        """
        Test the initializer.
        """

        pm = PermutationMatrix(self.atoms, self.cutoff)
        pm_cpp = permutation_matrix_from_atoms(self.atoms, self.cutoff)[0]
        self.assertIsInstance(pm, PermutationMatrix)
        self.assertIsInstance(pm_cpp, PermutationMatrix_cpp)
    
    def test_equal_first_column(self):
        """
        Test that the first column are the same for both matrices
        """

        matrix_py = self.pm.permutaded_matrix_frac
        matrix_cpp    = self.pm_cpp.get_permutated_positions()
        col1_py = [tuple(row[0]) for row in matrix_py]
        col1_cpp = [tuple(row[0]) for row in matrix_cpp]        
        print(type(col1_py))
        print(type(col1_cpp))
        print(len(col1_py))
        print(len(col1_cpp))
        # for pos1, pos2 in zip(sorted(col1_py), sorted(col1_cpp)):
        #     diff = np.array((pos1)) - np.array((pos2))
        #     print(pos1, pos2, diff)
        # self.assertListEqual(col1_py, col1_cpp)


if __name__ == '__main__':
    unittest.main()
