import unittest
from ase.build import bulk
from icet import Structure
from icet.core.lattice_site import LatticeSite
import random
import numpy as np


class TestStructure(unittest.TestCase):
    '''
    Container for test of the module functionality.
    Test is done againts ASE Atoms object
    '''

    def __init__(self, *args, **kwargs):
        super(TestStructure, self).__init__(*args, **kwargs)
        self.atoms = bulk("Al")
        self.noise = 1e-6

    def setUp(self):
        '''
        SetUp
        '''
        self.structure = Structure.from_atoms(self.atoms)

    def test_structure_from_atoms(self):
        '''
        Test that Structure object is returned form an
        ASE Atoms object
        '''
        pass

    def test_find_lattice_site_by_position_simple(self):
        """
        Tests finding lattice site by position, simple version using
        only one atom cell.

        1. Create a bunch of lattice sites all with index 0 and
        integer unitcell offsets
        2. convert these to x,y,z positions. Nothing strange so far
        3. Find lattice site from the position and assert that it should
           be equivalent to the original lattice site.
        """

        lattice_sites = []
        noise_position = []
        unit_cell_range = 1000
        for j in range(5000):
            offset = [random.randint(-unit_cell_range, unit_cell_range)
                      for i in range(3)]
            noise_position.append(
                [self.noise * random.uniform(-1, 1) for i in range(3)])
            lattice_sites.append(LatticeSite(0, offset))

        positions = []
        for i, site in enumerate(lattice_sites):
            # Get position with a little noise
            pos = self.structure.get_position(site)
            pos = pos + np.array(noise_position[i])
            positions.append(pos)
        for site, pos in zip(lattice_sites, positions):
            found_site = self.structure.find_lattice_site_by_position(pos)
            self.assertEqual(site, found_site)

    def test_find_lattice_site_by_position_medium(self):
        """
        Tests finding lattice site by position, medium version
        tests against hcp and user more than one atom in the basis
        1. Create a bunch of lattice sites all with index 0 and
        integer unitcell offsets
        2. convert these to x,y,z positions. Nothing strange so far
        3. Find lattice site from the position and assert that it should
           be equivalent to the original lattice site.
        """
        atoms = bulk("Au", 'hcp', a=2.0).repeat([3, 2, 5])
        structure = Structure.from_atoms(atoms)
        lattice_sites = []
        unit_cell_range = 1000
        noise_position = []

        for j in range(5000):
            offset = [random.randint(-unit_cell_range, unit_cell_range)
                      for i in range(3)]
            index = random.randint(0, len(atoms) - 1)
            noise_position.append(
                [self.noise * random.uniform(-1, 1) for i in range(3)])
            lattice_sites.append(LatticeSite(index, offset))

        positions = []
        for i, site in enumerate(lattice_sites):
            pos = structure.get_position(site)
            pos = pos + np.array(noise_position[i])
            positions.append(pos)
        for site, pos in zip(lattice_sites, positions):
            found_site = structure.find_lattice_site_by_position(pos)

            self.assertEqual(site, found_site)

    def test_find_lattice_site_by_position_hard(self):
        """
        Tests finding lattice site by position, hard version tests against hcp,
        many atoms in the basis AND pbc = [True, True, False] !
        1. Create a bunch of lattice sites all with index 0 and
        integer unitcell offsets
        2. convert these to x,y,z positions. Nothing strange so far
        3. Find lattice site from the position and assert that it should
           be equivalent to the original lattice site.
        """

        atoms = bulk("Au", 'hcp', a=2.0).repeat([3, 5, 5])

        # Set pbc false in Z-direction and add vacuum
        atoms.pbc = [True, True, False]
        atoms.center(30, axis=[2])
        structure = Structure.from_atoms(atoms)
        noise_position = []

        lattice_sites = []
        unit_cell_range = 100
        for j in range(500):
            offset = [random.randint(-unit_cell_range, unit_cell_range)
                      for i in range(3)]
            offset[2] = 0
            index = random.randint(0, len(atoms) - 1)
            noise_position.append(
                [self.noise * random.uniform(-1, 1) for i in range(3)])

            lattice_sites.append(LatticeSite(index, offset))

        positions = []
        for i, site in enumerate(lattice_sites):
            pos = structure.get_position(site)
            pos += np.array(noise_position[i])
            positions.append(pos)
        for site, pos in zip(lattice_sites, positions):
            found_site = structure.find_lattice_site_by_position(pos)
            self.assertEqual(site, found_site)

    def test_structure_to_atoms(self):
        '''
        Test that Structure is returned as an ASE
        Atoms object
        '''
        pass

    def test_repr_function(self):
        '''
        Test representation of a Structure object
        '''
        pass


if __name__ == '__main__':
    unittest.main()
