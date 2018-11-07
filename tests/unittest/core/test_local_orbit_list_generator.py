import unittest
from ase import Atoms
import numpy as np

from icet import OrbitList, Structure
from icet.core.local_orbit_list_generator import LocalOrbitListGenerator


class TestLocalOrbitListGenerator(unittest.TestCase):
    """
    Container for test of class functionality.
    """
    def __init__(self, *args, **kwargs):
        super(TestLocalOrbitListGenerator, self).__init__(*args, **kwargs)
        # corner case: tilted structure
        atoms = Atoms('Al',
                      positions=[[0., 0., 0.]],
                      cell=np.array([[0., 4., 4.],
                                     [4., 0., 4.],
                                     [4., 4., 0.]]),
                      pbc=[1, 1, 1])
        cutoffs = [4.2, 4.2]
        self.orbit_list = OrbitList(atoms, cutoffs)
        self.primitive = Structure.from_atoms(atoms)
        self.supercell = Structure.from_atoms(atoms.repeat(2))

    def shortDescription(self):
        """Silences unittest from printing the docstrings in test cases."""
        return None

    def setUp(self):
        """Instantiate class for each test case."""
        self.lol_gen = LocalOrbitListGenerator(self.orbit_list, self.supercell)

    def test_generate_local_orbit_list_from_index(self):
        """
        Tests that function generates an orbit list from
        an index of a specific offset of the supercell.
        """
        unique_offsets = self.lol_gen.get_unique_primcell_offsets()

        for index, offset in enumerate(unique_offsets):
            local_orbit_list = self.lol_gen.generate_local_orbit_list(index)
            for orbit_prim, orbit_super in zip(self.orbit_list.orbits,
                                               local_orbit_list.orbits):
                for sites_p, sites_s in zip(orbit_prim.representative_sites,
                                            orbit_super.representative_sites):
                    sites_p.unitcell_offset += offset
                    # pos_super = self.supercell.get_position(sites_s)
                    # pos_prim = self.primitive.get_position(sites_p)
                    print(sites_p, '|', sites_s)
                print('-------')

    def test_generate_local_orbit_list_from_offset(self):
        """
        Tests that function generates an orbit list by passing
        an index that corresponds to a specific offset of the supercell.
        """
        unique_offsets = self.lol_gen.get_unique_primcell_offsets()
        for offset in unique_offsets:
            local_orbit_list = self.lol_gen.generate_local_orbit_list(offset)
            for orbit_prim, orbit_super in zip(self.orbit_list.orbits,
                                               local_orbit_list.orbits):
                for sites_p, sites_s in zip(orbit_prim.representative_sites,
                                            orbit_super.representative_sites):
                    sites_p.unitcell_offset += offset
                    pos_super = self.supercell.get_position(sites_s)
                    pos_prim = self.primitive.get_position(sites_p)
                    print(pos_prim, '|', pos_super)

    def test_generate_full_orbit_list(self):
        """Tests that local orbit list of all unique offsets are returned."""
        full_orbit_list = self.lol_gen.generate_full_orbit_list()
        print(full_orbit_list)

    def test_clear(self):
        """
        Tests vector of unique offsets and primitive to supercell mapping
        are cleared.
        """
        self.lol_gen.generate_local_orbit_list(0)
        self.lol_gen.clear()
        offsets_count = self.lol_gen.get_unique_offsets_count()
        self.assertEqual(offsets_count, 0)
        mapping = self.lol_gen.get_primitive_to_supercell_map()
        self.assertEqual(len(mapping), 0)

    def test_unique_offset_count(self):
        """
        Tests number of unique offsets corresponds to the number of atoms
        in the supercell given that there is one atom in the primitive cell.
        """
        self.assertEqual(self.lol_gen.get_unique_offsets_count(),
                         len(self.supercell))

    def test_get_primitive_to_supercell_map(self):
        """
        Tests primitive to supercell mapping.
        """
        unique_offsets = self.lol_gen.get_unique_primcell_offsets()
        unique_offsets = [unique_offsets[-1]]

        for offset in unique_offsets:
            print('new offset: ', offset)
            self.lol_gen.generate_local_orbit_list(offset)
            mapping = self.lol_gen.get_primitive_to_supercell_map()
            for sites_prim, sites_super in mapping.items():
                pos_super = self.supercell.get_position(sites_super)
                pos_prim = self.primitive.get_position(sites_prim)
                print(sites_prim, '|', sites_super)
                print(pos_prim, '|', pos_super)

    def test_unique_primcell_offsets(self):
        """Tests uniqueness of the offsets of primitive cell."""
        offsets = self.lol_gen.get_unique_primcell_offsets()
        for i, offset in enumerate(offsets):
            for k in range(i+1, len(offsets)):
                self.assertFalse(np.all(np.isclose(offset, offsets[k])))


if __name__ == '__main__':
    unittest.main()
