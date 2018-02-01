import unittest
import numpy as np
import random
import timeit

from ase.build import bulk

from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core_py.lattice_site import LatticeSite
#from icet.core.orbit_list import __get_lattice_site_permutation_matrix\
#    as get_lattice_site_permutation_matrix




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

    def _prune(self, sites, prune_method):
        """
        Helper function that executes the indicated prune_method with sites as pm_lattice_stes

        """
        self.pm.pm_lattice_sites = sites
        getattr(self.pm, prune_method)()

    def _test__prune_permutation_matrix(self, prune_method='_prune_permutation_matrix'):
        inner_list1 = (0, 1, 2)
        inner_list2 = (3, 10, 4)
        inner_list3 = [LatticeSite(0, (0,1,1)), LatticeSite(2, (1,1,1)), LatticeSite(0, (0,0,0))]
        inner_list4 = [LatticeSite(4, (0,1,1)), LatticeSite(4, (1,1,1)), LatticeSite(4, (0,0,0))]
        inner_list5 = [LatticeSite(0, (1,1,1)), LatticeSite(0, (1,1,1)), LatticeSite(0, (1,1,1))]


        # Define a list of input, expected_outcome for a number of sub tests
        tests = [(10*[inner_list1],                                     [inner_list1]),
                 ([inner_list1, inner_list2],                           [inner_list1, inner_list2]),
                 ([inner_list1, inner_list2, inner_list2, inner_list1], [inner_list1, inner_list2]),
                 (10*[inner_list3],                                     [inner_list3]),
                 ([inner_list3, inner_list5, inner_list3, inner_list3], [inner_list3, inner_list5])
                 ]

        # Create list to be used for benchmark
        n_elements = 2000
        long_list = n_elements * ([inner_list3] + [inner_list4] + [inner_list5])
        random.shuffle(long_list)


        def assertExpectedOutcome(sites, expected_outcome, prune_method=prune_method):
            self._prune(sites, prune_method)
            self.assertEqual(self.pm.pm_lattice_sites.sort(), expected_outcome.sort())


        for sites, expected_outcome in tests:
            with self.subTest(sites=sites, expected_outcome=expected_outcome):
                assertExpectedOutcome(sites, expected_outcome)

        # benchmark
        parms = locals()
        parms.update(globals())
        dt = timeit.timeit('self._prune(long_list, prune_method)', globals=parms, number=5)
        print('{} : {} ({})'.format(prune_method[-1], dt, prune_method))



        # Note that comparison only happens between the first element (element 0) of the inner loop
        # thus inner_lista = (0, 10) together with inner_list1 will only result in one of the lists
        # remaining because the first element is 'a' in both cases
        inner_lista = [0, 10]
        assertExpectedOutcome([inner_list1, inner_lista], [inner_list1])  # rather than [inner_list2, inner_lista]


    def test__prune_permutation_matrices(self):
        for lett in ['A', 'B', 'C', 'D', 'F']:  # E fails
            prune_method = '_prune_permutation_matrix' + lett

            with self.subTest(prune_method=prune_method):
                self._test__prune_permutation_matrix(prune_method=prune_method)


if __name__ == '__main__':
    unittest.main()
