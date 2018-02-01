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

    def setUp(self):
        self.atoms = bulk("Al")
        self.cutoff = 5.0

    def test__init__(self):
        """
        Test that we can initialize the PermutationMatrix
        """
        try:
            pm = PermutationMatrix(self.atoms, self.cutoff)
        except Exception as inst:
            self.fail('Failed with exception {}'.format(inst))



if __name__ == '__main__':
    unittest.main()
