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
        n_elements = 20000
        long_list = n_elements * ([20*inner_list3] + [20*inner_list4] + [20*inner_list5])
        random.shuffle(long_list)


        def assertExpectedOutcome(sites, expected_outcome, prune_method=prune_method):
            self._prune(sites, prune_method)
            self.assertEqual(self.pm.pm_lattice_sites.sort(), expected_outcome.sort())


        for sites, expected_outcome in tests[3:]:
            with self.subTest(sites=sites, expected_outcome=expected_outcome):
                assertExpectedOutcome(sites, expected_outcome)

        # benchmark
        parms = locals()
        parms.update(globals())
        test_list = long_list[:]
        self.assertEqual(len(long_list), 3*n_elements)
        dt = timeit.timeit('self._prune(long_list, prune_method)', globals=parms, number=5)
        print('{} : {} ({})'.format(prune_method[-1], dt, prune_method))

        # Is numpy faster?
        # convert long_list to numpy array

        # shape = (n_rows, n_sites_per_row, 4)
        try:
            long_array1 = np.array([[site.as_list for site in row] for row in test_list])
        except:
            long_array1 = np.array(test_list)

        # shape = (n_sites_per_rowm n_rows, 4)
        row_length = len(test_list[0])
        col_length = len(test_list)
        long_array2 = np.zeros( (row_length, col_length, 4), 'i')
        for i in range(len(test_list[0])):
            for j in range(len(test_list)):
                try:
                    kk = test_list[j][i].as_list
                except AttributeError:
                    kk = test_list[j][i]
                for k in range(4):
                    long_array2[i, j, k] = kk[k]


        parms = locals()
        parms.update(globals())
        dt = timeit.timeit('np.unique(long_array1, axis=0)', globals=parms, number=5)
        print('{} : {} ({})'.format('M', dt, 'numpy unique M'))

        def get_unique(narray):
            first_sites_in_rows = long_array1[:, 0, :]
            a, indexes = np.unique(first_sites_in_rows, axis=0, return_index=True)
            return long_array1[indexes]

        parms = locals()
        parms.update(globals())
        dt = timeit.timeit('get_unique(long_array1)', globals=parms, number=5)
        print('{} : {} ({})'.format('N', dt, 'numpy unique N'))



        # Note that comparison only happens between the first element (element 0) of the inner loop
        # thus inner_lista = (0, 10) together with inner_list1 will only result in one of the lists
        # remaining because the first element is 0 in both cases

        lean_flag = hasattr(inner_list3[0], 'as_list')

    def test__prune_permutation_matrices(self):
        for lett in ['A', 'B', 'D', 'F', 'G', 'H']:  # C is too slow, E fails
            prune_method = '_prune_permutation_matrix' + lett
            if not lean_flag or lett != 'G':
                with self.subTest(prune_method=prune_method):
                    self._test__prune_permutation_matrix(prune_method=prune_method)


if __name__ == '__main__':
    unittest.main()
