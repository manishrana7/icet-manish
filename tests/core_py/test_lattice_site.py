from icet.core_py.lattice_site import LatticeSite as LatticeSite


import itertools
import unittest


class TestLatticeSite(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        super(TestLatticeSite, self).__init__(*args, **kwargs)

        self.indices = [i for i in range(9)]
        self.offsets = []
        cartesian_product_lists = [[0, 1], [0, 1], [0, 1]]
        for element in itertools.product(*cartesian_product_lists):
            print(element)
            self.offsets.append(element)

    def setUp(self):
        """
        Setup
        """
        self.lattice_sites = []
        for index, offset in zip(self.indices, self.offsets):
            lattice_site = LatticeSite(index, offset)
            self.lattice_sites.append(lattice_site)

    def test_printing(self):
        """
        Test printing a LatticeSite
        """
        print(self.lattice_sites[0])


if __name__ == '__main__':
    unittest.main()
