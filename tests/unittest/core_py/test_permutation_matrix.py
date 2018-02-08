import unittest
import numpy as np
import random
import timeit

from ase.build import bulk

from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core_py.lattice_site import LatticeSite

# from icet.core_py.lattice_site import LatticeSiteLean as LatticeSite

#from icet.core.orbit_list import __get_lattice_site_permutation_matrix\
#    as get_lattice_site_permutation_matrix


lean_flag = False

if lean_flag:
    from icet.core_py.lattice_site import LatticeSiteLean as LatticeSite
else:
    from icet.core_py.lattice_site import LatticeSite


class TestPermutationMatrix(unittest.TestCase):
    """
    Test the python implementation of permutatation matrix
    """

    def __init__(self, *args, **kwargs):
        super(TestPermutationMatrix, self).__init__(*args, **kwargs)
        self.al_atoms = bulk("Al")
        self.cutoff = 3.0


    def setUp(self):
        self.pm = PermutationMatrix(self.al_atoms, self.cutoff)
        self.lattice_sites = self.pm.pm_lattice_sites


    def test_type_of_pm_lattice_sites(self):
        """
        Test that the permutation matrix is a list
        """
        self.assertTrue(isinstance(self.lattice_sites, list))
        self.assertEqual(self.lattice_sites.__class__.__name__, 'list')


    def test_len_of_pm_lattice_sites(self):
        """
        Test that the permutation matrix has the correct length
        """
        self.assertEqual(13, len(self.lattice_sites))


    def test_type_of_pm_lattice_sites_row(self):
        """
        Test that a is of type list
        """
        self.assertTrue(isinstance(self.lattice_sites[0], list))
        self.assertEqual(self.lattice_sites[0].__class__.__name__, 'list')


    def test_len_of_pm_lattice_sites_row(self):
        """
        Test that a row has the correct length
        """
        self.assertEqual(48, len(self.lattice_sites[0]))


    def test_lattice_sites_element(self):
        """ Test that first element in the permutation matrix is correct"""
        site = self.lattice_sites[0][0]
        self.assertTrue(isinstance(site, LatticeSite))


    def test_lattice_sites_all_rows(self):
        """
        Test that all rows have the same length
        """
        n = len(self.lattice_sites[0])
        diff_lens = [row for row in self.lattice_sites if len(row) != n]
        self.assertEqual(diff_lens, [])


if __name__ == '__main__':
    unittest.main()
