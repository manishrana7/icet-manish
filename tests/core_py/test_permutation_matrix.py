from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core.permutation_map import permutation_matrix_from_atoms
from icet.core.permutation_map import PermutationMap as PermutationMatrix_cpp
from ase.build import bulk
from icet.core.orbit_list import __get_lattice_site_permutation_matrix\
    as get_lattice_site_permutation_matrix

import unittest

import numpy as np


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
        matrix_cpp = sorted(self.pm_lattice_sites_cpp)
        col1_py = [row[0] for row in matrix_py]
        col1_cpp = [row[0] for row in matrix_cpp]
        for lattice_site_py, lattice_site_cpp in zip((col1_py),
                                                     (col1_cpp)):
            self.assertEqual(lattice_site_py, lattice_site_cpp)
        self.assertListEqual((col1_py), (col1_cpp))

    def test_all_matrix_elements_equal(self):
        """
        Test that all elemetns of both Lattice site permutation
        matrices are identical.
        """
        pm_cpp = sorted(self.pm_lattice_sites_cpp)
        for row_py, row_cpp in zip(self.pm.pm_lattice_sites, pm_cpp):
            self.assertEqual(len(row_py), len(row_cpp))
            for element_py, element_cpp in zip(row_py, row_cpp):
                self.assertEqual(element_py, element_cpp)

    def test_primitive_structure(self):
        """
        Test primitive structure property by comparing to
        the other versions primitive structure.
        """
        for position_py, position_cpp in zip(
                self.pm.primitive_structure.positions,
                self.prim_structure_cpp.positions):
            self.assertAlmostEqual(np.linalg.norm(
                position_cpp - position_py), 0)

        self.assertEqual(self.pm.primitive_structure.cell.all(),
                         self.prim_structure_cpp.cell.all())


if __name__ == '__main__':
    unittest.main()
