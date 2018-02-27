import unittest

from ase.build import bulk
from icet import Structure
from icet.core.lattice_site import LatticeSite
from icet.core.neighbor_list import NeighborList
from icet.core.many_body_neighbor_list import (
    ManyBodyNeighborList, get_all_lattice_neighbors)


class TestManyBodyNeighborList(unittest.TestCase):
    """
    Container for test of the module functionality.
    """

    def __init__(self, *args, **kwargs):
        super(TestManyBodyNeighborList, self).__init__(*args, **kwargs)

        self.atoms = bulk('Al').repeat(2)
        # self.atoms = bulk('Ni', 'hcp', a=3.0).repeat([2, 2, 1])
        self.cutoffs = [5.0, 5.0]

    def setUp(self):
        """
        SetUp
        """
        self.mbnl = ManyBodyNeighborList()

        self.neighbor_lists = []
        for cutoff in self.cutoffs:
            structure = Structure.from_atoms(self.atoms)
            nl = NeighborList(cutoff)
            nl.build(structure)
            self.neighbor_lists.append(nl)

    def test_build(self):
        """
        Test build.
        """
        for index in range(len(self.atoms)):
            self.mbnl.build(self.neighbor_lists, index, True)

    def test_bothways_true(self):
        """
        Build the mbnl with bothways = True and assert that
        each index in the atoms object have the same number
        of neighbors.
        """
        mbnl_size = len(self.mbnl.build(self.neighbor_lists, 0, True))
        for index in range(1, len(self.atoms)):
            self.assertEqual(mbnl_size, len(self.mbnl.build(
                self.neighbor_lists, index, True)))

    def test_bothways_false(self):
        """
        Build the mbnl with bothways = False and assert that
        mbnl built on the first index in the atoms object do not
        have the same number of neighbors as the other atoms.
        """
        mbnl_size = len(self.mbnl.build(self.neighbor_lists, 0, False))
        for index in range(1, len(self.atoms)):
            self.assertNotEqual(mbnl_size, len(self.mbnl.build(
                self.neighbor_lists, index, False)))

    def test_singlets(self):
        """
        Test singlet.
        """
        for index in range(len(self.atoms)):
            target = tuple(([LatticeSite(index, [0., 0., 0.])], []))
            neighbors = self.mbnl.build(self.neighbor_lists, index, False)
            self.assertEqual(neighbors[0], target)

    def test_pairs(self):
        """
        Test pairs.
        """
        index = 0
        nl_neighbors = self.neighbor_lists[0].get_neighbors(0)
        target = tuple(([LatticeSite(index, [0., 0., 0.])], nl_neighbors))
        neighbors = self.mbnl.build(self.neighbor_lists, index, False)
        self.assertEqual(neighbors[1], target)

    def test_calculate_intersections(self):
        """
        Test intersection between two list of neighbors.
        """
        lattice_sites = []
        lattice_sites.append(LatticeSite(0, [0, 0, 0]))
        lattice_sites.append(LatticeSite(0, [1, 0, 0]))
        lattice_sites.append(LatticeSite(1, [0, 0, 0]))
        lattice_sites.append(LatticeSite(3, [0, 0, 0]))

        lattice_sites2 = []
        lattice_sites2.append(LatticeSite(0, [0, 0, 0]))
        lattice_sites2.append(LatticeSite(0, [1, 0, 0]))

        intersection = self.mbnl.calculate_intersection(
            lattice_sites, lattice_sites2)

        self.assertEqual(sorted(intersection), [LatticeSite(
            0, [0, 0, 0]), LatticeSite(0, [1, 0, 0])])

    def test_get_all_lattice_neighbors(self):
        """
        Test get_all_lattice_neighbors functionality.
        """
        neighbors = get_all_lattice_neighbors(self.atoms,
                                              self.neighbor_lists,
                                              self.cutoffs)
        neighbors2 = get_all_lattice_neighbors(atoms=self.atoms,
                                               cutoffs=self.cutoffs)
        for neighbor, neighbor2 in zip(neighbors, neighbors2):
            self.assertEqual(neighbor, neighbor2)

    def test_mbnl_non_pbc(self):
        """
        Test many-body neighbor list for non-pbc structure gives
        less neighbors under the same cutoffs.
        """
        atoms = self.atoms.copy()
        atoms.set_pbc([False])
        neighbor_lists = []
        for cutoff in self.cutoffs:
            structure = Structure.from_atoms(atoms)
            nl = NeighborList(cutoff)
            nl.build(structure)
            neighbor_lists.append(nl)

        neighbors_non_pbc = get_all_lattice_neighbors(atoms,
                                                      neighbor_lists,
                                                      self.cutoffs)

        neighbors = get_all_lattice_neighbors(self.atoms,
                                              self.neighbor_lists,
                                              self.cutoffs)

        self.assertLess(len(neighbors_non_pbc), len(neighbors))


if __name__ == '__main__':
    unittest.main()
