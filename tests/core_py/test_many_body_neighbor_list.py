import unittest

from ase.build import bulk

from icet.core_py.many_body_neighbor_list import ManyBodyNeighborList

class TestManyBodyNeighborList(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestManyBodyNeighborList, self).__init__(*args, **kwargs)
        
        self.atoms_prim = bulk("Al")
        self.atoms = bulk("Al").repeat(2)
        self.cutoffs = [5, 5]

    def setUp(self):
        """
        Setup before each test.        
        """
        self.mbnl = ManyBodyNeighborList()

    def test_build(self):
        """
        Test that a simple build works
        """
        index = 0
        self.mbnl.build(self.neighbor_lists, index, bothways=False)


    



if __name__ == '__main__':
    unittest.main()
