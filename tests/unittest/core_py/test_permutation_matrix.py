import unittest
import numpy as np

from ase.build import bulk

from icet.core_py.permutation_matrix import PermutationMatrix
#from icet.core.orbit_list import __get_lattice_site_permutation_matrix\
#    as get_lattice_site_permutation_matrix




class TestPermutationMatrix(unittest.TestCase):
    """
    Test the python implementation of permutatation matrix
    """

    def __init__(self, *args, **kwargs):
        super(TestPermutationMatrix, self).__init__(*args, **kwargs)
        self.atoms = bulk("Al")
        self.cutoff = 5.0

    def test__init__(self, is_forced_to_fail=True):
        """
        Test that we can instantiate the PermutationMatrix
        """
        pm = PermutationMatrix(self.atoms, self.cutoff)


    def test_shape_of_pm_lattice_sites(self):
        """
        Test that the 'pm_lattice_sites' attribute to PermutationMatrix has the right shape and dimensions
        """
        pm = PermutationMatrix(self.atoms, self.cutoff)
        lattice_sites = pm.pm_lattice_sites



if __name__ == '__main__':
    unittest.main()
