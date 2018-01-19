import unittest

from ase.build import bulk
from ase.neighborlist import NeighborList
from icet.core.neighbor_list import NeighborList as NeighborList_cpp
from icet.core_py.many_body_neighbor_list import ManyBodyNeighborList
from icet import Structure
from icet.core.many_body_neighbor_list import (
    ManyBodyNeighborList as ManyBodyNeighborList_cpp)


class TestManyBodyNeighborList(unittest.TestCase):
    """
        TODO
        ----
        * compare each latticesite or list of lattice site in each test
          since right now you might not test true equivalence between versions.
    """

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
        self.mbnl_cpp = ManyBodyNeighborList_cpp()
        self.neighbor_lists = []
        self.neighbor_lists_cpp = []
        for co in self.cutoffs:
            ase_nl = NeighborList(len(self.atoms) * [co / 2], skin=1e-8,
                                  bothways=True, self_interaction=False)

            ase_nl.update(self.atoms)

            structure = Structure.from_atoms(self.atoms)
            nl = NeighborList_cpp(co)
            nl.build(structure)

            self.neighbor_lists.append(ase_nl)
            self.neighbor_lists_cpp.append(nl)

    def test_build(self):
        """
        Test that a simple build works
        """
        for index in range(len(self.atoms)):
            mbnl_py = self.mbnl.build(
                self.neighbor_lists, index, bothways=False)
            mbnl_cpp = self.mbnl_cpp.build(
                self.neighbor_lists_cpp, index, False)
            self.assertEqual(len(mbnl_py), len(mbnl_cpp))

    def test_bothways_true(self):
        """
        Build the mbnl with bothways = True and
        assert that each index in the atoms object
        have the same number of neighbors


        """

        mbnl_size = len(self.mbnl.build(self.neighbor_lists, 0, bothways=True))
        mbnl_size_cpp = len(self.mbnl_cpp.build(
            self.neighbor_lists_cpp, 0, True))

        for index in range(len(self.atoms)):
            self.assertEqual(mbnl_size, len(self.mbnl.build(
                self.neighbor_lists, index, bothways=True)))
            self.assertEqual(mbnl_size_cpp, len(self.mbnl_cpp.build(
                self.neighbor_lists_cpp, index,
                True)))

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
        mbnl_size_cpp = len(self.mbnl_cpp.build(
            self.neighbor_lists_cpp, 0, False))
        for index in range(1, len(self.atoms)):
            self.assertNotEqual(mbnl_size, len(self.mbnl.build(
                                self.neighbor_lists, index,
                                bothways=False)))
            self.assertNotEqual(mbnl_size_cpp, len(self.mbnl_cpp.build(
                                self.neighbor_lists_cpp, index,
                                False)))

        # compare to cpp
        for index in range(len(self.atoms)):
            mbnl_py = self.mbnl.build(
                self.neighbor_lists, index, bothways=False)
            mbnl_cpp = self.mbnl_cpp.build(
                self.neighbor_lists_cpp, index, False)
            for lat_site_cpp, lat_site in zip(mbnl_py, mbnl_cpp):
                self.assertEqual(lat_site_cpp[0], lat_site[0])
                self.assertEqual(lat_site_cpp[1], lat_site[1])


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
