import unittest
from ase.build import bulk
from ase import Atoms
import numpy as np

from icet.core.structure import Structure
from icet.core.orbit_list import create_orbit_list
from icet.core.local_orbit_list_generator import LocalOrbitListGenerator

class TestLocalOrbitListGenerator(unittest.TestCase):
    """
    Container for test of class functionality.
    """

    def __init__(self, *args, **kwargs):
        super(TestLocalOrbitListGenerator, self).__init__(*args, **kwargs)
        atoms = Atoms('Al',
                      positions=[[0., 0., 0.]],
                      cell= np.array([[0., 4., 1.],
                                      [4., 0., 1.],
                                      [4., 4., 0.]]), 
                      pbc=[1, 1, 1])
        cutoffs = [3.2, 3.2]
        # Build an orbit list and supercell
        self.orbit_list = create_orbit_list(atoms, cutoffs)
        self.structure = Structure.from_atoms(atoms.repeat(2))

    def setUp(self):
        """
        Instantiate class for each test case.
        """
        self.lolg = LocalOrbitListGenerator(self.orbit_list, self.structure)

    def test_generate_local_orbit_list(self):
        """
        Tests that function generates an orbit list by passing
        an index that corresponds to a specific offset of the supercell.
        """
        index = 0
        local_orbit_list = self.lolg.generate_local_orbit_list(index)
        for orbit in local_orbit_list.get_orbit_list():
            print(orbit.get_representative_sites().__str__())

    @unittest.SkipTest
    def test_generate_full_orbit_list(self):
        """
        Tests functionality.
        """
        full_orbit_list = self.lolg.generate_full_orbit_list()
        for orbit in full_orbit_list.get_orbit_list():
            print(orbit.get_representative_sites().__str__())

    def test_clear(self):
        """
        Tests whatever is cleared by calling this function.
        """
        self.lolg.clear()
        index=0
        with self.assertRaises(IndexError):
            self.lolg.generate_local_orbit_list(index)

    def test_unique_offset_count(self):
        """
        Tests that unique offsets obtained by mapping supercell sites
        to primitive cell.
        """
        unique_offsets =  self.lolg.get_unique_offsets_count()
        print(unique_offsets)

    def test_get_primitive_to_supercell_map(self):
        """
        Tests functionality.
        """
        mapping = self.lolg.get_primitive_to_supercell_map()
        print(mapping)

    def test_unique_primcell_offsets(self):
        """
        Tests unique offset of primitive cell.
        """
        unique_offsets = self.lolg.get_unique_primcell_offsets()
        print(unique_offsets)


if __name__ == '__main__':
    unittest.main()
