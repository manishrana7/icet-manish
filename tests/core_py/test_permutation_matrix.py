from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core.permutation_map import permutation_matrix_from_atoms
from icet.core.permutation_map import PermutationMap as PermutationMatrix_cpp
from ase.build import bulk
from icet.core.orbit_list import __get_lattice_site_permutation_matrix\
    as get_lattice_site_permutation_matrix

import unittest


class TestPermutationMatrix(unittest.TestCase):
    '''
    Test the python implementation of permutation matrix
    '''

    def __init__(self, *args, **kwargs):
        super(TestPermutationMatrix, self).__init__(*args, **kwargs)

        self.atoms = bulk("Al")
        self.cutoff = 5.0

    def setUp(self):
        """
        Setup.
        """
        self.pm = PermutationMatrix(self.atoms, self.cutoff)
        self.pm_cpp, self.prim_structure_cpp, _ = \
            permutation_matrix_from_atoms(
                self.atoms, self.cutoff)
        self.pm_lattice_sites_cpp = \
            get_lattice_site_permutation_matrix(self.prim_structure_cpp,
                                                self.pm_cpp,
                                                prune=True)

    def test_init(self):
        """
        Test the initializer.
        """

        pm = PermutationMatrix(self.atoms, self.cutoff)
        pm_cpp, prim_structure_cpp, _ = permutation_matrix_from_atoms(
            self.atoms, self.cutoff)
        self.assertIsInstance(pm, PermutationMatrix)
        self.assertIsInstance(pm_cpp, PermutationMatrix_cpp)

    def test_equal_first_column(self):
        """
        Test that the first column are the same for both matrices
        """

        matrix_py = self.pm.pm_lattice_sites
        matrix_cpp = self.pm_lattice_sites_cpp
        col1_py = [row[0] for row in matrix_py]
        col1_cpp = [row[0] for row in matrix_cpp]
        for lattice_site_py, lattice_site_cpp in zip(sorted(col1_py),
                                                     sorted(col1_cpp)):
            self.assertEqual(lattice_site_py, lattice_site_cpp)
        self.assertListEqual(sorted(col1_py), sorted(col1_cpp))


if __name__ == '__main__':
    unittest.main()
