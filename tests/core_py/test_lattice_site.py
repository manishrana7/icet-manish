from icet.core.lattice_site import LatticeSite as LatticeSite_cpp
from icet.core_py.lattice_site import LatticeSite as LatticeSite


import itertools
import unittest


class TestLatticeSite(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        super(TestLatticeSite, self).__init__(*args, **kwargs)

        self.indices = [i for i in range(8)]
        self.offsets = []
        cartesian_product_lists = [[0, 1], [0, 1], [0, 1]]
        for element in itertools.product(*cartesian_product_lists):
            self.offsets.append(list(element))

    def setUp(self):
        """
        Setup
        """
        self.lattice_sites = []
        self.lattice_sites_cpp = []
        for index, offset in zip(self.indices, self.offsets):
            lattice_site = LatticeSite(index, offset)
            lattice_site_cpp = LatticeSite_cpp(index, offset)
            self.lattice_sites.append(lattice_site)
            self.lattice_sites_cpp.append(lattice_site_cpp)

    def test_printing(self):
        """
        Test printing a LatticeSite
        """
        print(self.lattice_sites[0])
        print(self.lattice_sites_cpp[0])
        
        self.assertEqual(self.lattice_sites[0].__str__(),self.lattice_sites_cpp[0].__str__())

    def test_sorting(self):
        """
        Test sorting the lattice sites
        """
        self.lattice_sites.sort()
        self.assertEqual(self.lattice_sites[0].index, 0)
        self.assertEqual(self.lattice_sites[0].offset, [0, 0, 0])

        self.lattice_sites.sort(reverse=True)
        self.assertEqual(self.lattice_sites[0].index, 7)
        self.assertEqual(self.lattice_sites[0].offset, [1, 1, 1])

    def test_lt(self):
        """
        Test less than operator.
        """
        self.assertLess(self.lattice_sites[0], self.lattice_sites[1])

    def test_eq(self):
        """
        Test eq operator
        """
        index = 15245345345
        offset = [-234234, 32423423, 235567567]
        lattice_site = LatticeSite(index, offset)
        lattice_site_other = LatticeSite(index, offset)
        self.assertEqual(lattice_site, lattice_site_other)

    def test_hash(self):
        """
        Test hash function
        """

        index = 15245345345
        offset = [-234234, 32423423, 235567567]
        lattice_site = LatticeSite(index, offset)
        lattice_site_other = LatticeSite(index, offset)
        self.assertEqual(lattice_site.__hash__(),
                         lattice_site_other.__hash__())



if __name__ == '__main__':
    unittest.main()
