import unittest

from ase.neighborlist import NeighborList as ASENeighborList
from icet import Structure
from _icet import NeighborList
from icet.core.neighbor_list import get_neighbor_lists
from ase.build import bulk
import numpy as np


class TestNeighborList(unittest.TestCase):
    """
    Container for test of the module functionality.
    """
    def __init__(self, *args, **kwargs):
        super(TestNeighborList, self).__init__(*args, **kwargs)

        self.atoms = bulk("NaCl", 'rocksalt', a=1.0).repeat(2)
        self.cutoff = 1.4
        self.structure = Structure.from_atoms(self.atoms)
        self.ase_nl = ASENeighborList(len(self.atoms)*[self.cutoff/2],
                                      skin=1e-8, bothways=True,
                                      self_interaction=False)
        self.ase_nl.update(self.atoms)
        self.ase_indices, self.ase_offsets = self.ase_nl.get_neighbors(0)

    def setUp(self):
        """
        SetUp before each test.
        """
        self.nl = NeighborList(self.cutoff)
        self.nl.build(self.structure)
        self.neighbors = self.nl.get_neighbors(0)
        self.indices = []
        self.offsets = []
        for ngb in self.neighbors:
            self.indices.append(ngb.index)
            self.offsets.append(ngb.unitcell_offset)

    def test_build(self):
        """
        Test build gives the same number of neighbors
        as ASE neigborlist.
        """
        for index in range(len(self.atoms)):
            neighbors = self.nl.get_neighbors(index)
            ase_neighbors = self.ase_nl.get_neighbors(index)
            self.assertEqual(len(neighbors), len(ase_neighbors[0]))

    def test_indices(self):
        """
        Test indices are equal to those in ASE neighbor list.
        """
        self.assertEqual(len(self.indices), len(self.ase_indices))
        for index in self.indices:
            self.assertIn(index, self.ase_indices)

    def test_offsets(self):
        """
        Test offsets are equal to those in ASE neighbor list.
        """
        self.assertEqual(len(self.offsets), len(self.ase_offsets))
        for offset in self.offsets:
            self.assertIn(offset, self.ase_offsets)

    def test_equivalent_indices(self):
        """
        Test equivalent indices.
        """
        for index, offset in zip(self.indices, self.offsets):
            equiv_indices = [i for i, ase_offset in enumerate(self.ase_offsets)
                             if self.ase_indices[i] == index and
                             (ase_offset == offset).all()]
            self.assertEqual(len(equiv_indices), 1)
            self.assertEqual(self.ase_indices[equiv_indices[0]], index)

    def test_neighbors_lists(self):
        """
        Test neighbor lists againts ASE neighborlist.
        """
        for i in range(len(self.atoms)):
            index = [ngb.index for ngb in self.nl.get_neighbors(i)]
            offset = [ngb.unitcell_offset for ngb in self.nl.get_neighbors(i)]
            pos = self.atoms.positions[index] + np.dot(offset,
                                                       self.atoms.get_cell())
            index2, offset2 = self.ase_nl.get_neighbors(i)
            pos2 = self.atoms.positions[index2] + np.dot(offset2,
                                                         self.atoms.get_cell())
            self.assertCountEqual(pos2.tolist(), pos.tolist())

    def test_get_neighbors_lists(self):
        """
        Test get_neighbor_lists functionality.
        """
        nl = get_neighbor_lists(self.structure, [self.cutoff] * 4)
        self.assertEqual(len(nl), 4)
        self.assertEqual(len(nl[0]), len(self.nl))


if __name__ == '__main__':
    unittest.main()
