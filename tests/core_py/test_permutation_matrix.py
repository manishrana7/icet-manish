from icet.core_py.permutation_matrix import PermutationMatrix

from ase.build import bulk

import unittest


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

    def test_init(self):
        """
        Test the initializer.
        """

        pm = PermutationMatrix(self.atoms, self.cutoff)
        self.assertIsInstance(pm, PermutationMatrix)


if __name__ == '__main__':
    unittest.main()
