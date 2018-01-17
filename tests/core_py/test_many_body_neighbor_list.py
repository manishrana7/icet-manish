import unittest

from ase.build import bulk
from ase.neighborlist import NeighborList

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
        self.neighbor_lists = []
        for co in self.cutoffs:
            ase_nl = NeighborList(len(self.atoms) * [co / 2], skin=1e-8,
                                  bothways=True, self_interaction=False)
            ase_nl.update(self.atoms)
            self.neighbor_lists.append(ase_nl)

    def test_build(self):
        """
        Test that a simple build works
        """
        index = 0
        self.mbnl.build(self.neighbor_lists, index, bothways=False)

    def test_bothways_true(self):
        """
        Build the mbnl with bothways = True and
        assert that each index in the atoms object
        have the same number of neighbors
        """

        mbnl_size = len(self.mbnl.build(self.neighbor_lists, 0, bothways=True))
        for index in range(len(self.atoms)):
            self.assertEqual(mbnl_size, len(self.mbnl.build(
                self.neighbor_lists, index, bothways=True)))

    def test_bothways_false(self):
        """
        Build the mbnl with bothways = False and
        assert that mbnl built on the first 
        index in the atoms object do not
        have the same number of neighbors as 
        the other atoms.
        """

        mbnl_size = len(self.mbnl.build(
            self.neighbor_lists, 0, bothways=False))
        for index in range(1, len(self.atoms)):
            self.assertNotEqual(mbnl_size, len(self.mbnl.build(
                self.neighbor_lists, index, bothways=False)))


class BruteForceMBNL(object):
    """
    Builds an mbnl object from brute force.
    This is used to validate the more clever implemented
    mbnl.
    """
    def __init__():
        pass
    def build(self, neighbor_lists, index, bothways=False):
        """
        Will do something like this:
        mbnl = []
        for i in atoms:
            for j in neighbors[i]:
                for k in neighbors[j]:
                    if k in neigbhors[i]:
                    mbnl.append(i,j,k)
        """


if __name__ == '__main__':
    unittest.main()
