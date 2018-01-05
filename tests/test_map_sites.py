#!/usr/bin/env python3

import unittest

from icetdev.tools import map_structure_to_reference
from icetdev.tools.map_sites import (_get_scaled_cell,
                                     _get_transformation_matrix,
                                     _rescale_structures)
import numpy as np


class TestMapSites(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def __init__(self, *args, **kwargs):
        from ase import Atom
        from ase.build import bulk

        reference = bulk('Au', a=4.0)
        reference.append(Atom('H', (2, 2, 2)))
        self.reference = reference

        atoms = reference.repeat(2)
        for i in [0, 4, 6, 2, 7, 3]:
            if atoms[i].symbol == 'Au':
                atoms[i].symbol = 'Pd'
            elif atoms[i].symbol == 'H':
                del atoms[i]
        atoms.rattle(0.1)
        atoms.set_cell(atoms.cell * 1.01, scale_atoms=True)
        self.atoms = atoms

        super(TestMapSites, self).__init__(*args, **kwargs)

    def test_get_scaled_cell(self):
        '''
        Testing that the cell scaling retrieval works.
        '''
        modcell = _get_scaled_cell(self.atoms, self.reference,
                                   vacancy_type='V',
                                   inert_species=['Au', 'Pd'])
        target = np.array([[0.,  4.,  4.],
                           [4.,  0.,  4.],
                           [4.,  4.,  0.]])
        self.assertTrue(np.allclose(modcell, target))

        modcell = _get_scaled_cell(self.atoms, self.reference,
                                   vacancy_type='V')
        target = np.array([[0.,    4.04,  4.04],
                           [4.04,  0.,    4.04],
                           [4.04,  4.04,  0.]])
        self.assertTrue(np.allclose(modcell, target))

        modcell = _get_scaled_cell(self.atoms, self.reference)
        target = np.array([[0.,          3.82586237,  3.82586237],
                           [3.82586237,  0.,          3.82586237],
                           [3.82586237,  3.82586237,  0.]])
        self.assertTrue(np.allclose(modcell, target))

    def test_get_transformation_matrix(self):
        '''
        Testing that transformation matrix calculation works.
        '''
        P = _get_transformation_matrix(self.atoms.cell, self.reference.cell)
        target = np.array([[2.,  0.,  0.],
                           [0.,  2.,  0.],
                           [0.,  0.,  2.]])
        self.assertTrue(np.allclose(P, target))

    def test_rescale_structures(self):
        '''
        Testing that the structure rescaling works.
        '''
        P = np.array([[2.,  0.,  0.],
                      [0.,  2.,  0.],
                      [0.,  0.,  2.]])
        scaled_structure, ideal_supercell = \
            _rescale_structures(self.atoms, self.reference, P)
        self.assertEqual(len(scaled_structure), 14)
        self.assertEqual(len(ideal_supercell), 16)
        self.assertTrue(np.allclose(scaled_structure.cell,
                                    ideal_supercell.cell))

    def test_map_structure_to_reference(self):
        '''
        Testing that the mapping works. This function will also invoke the
        other ones.
        '''
        mapped, r_max, r_av = \
            map_structure_to_reference(self.atoms, self.reference,
                                       0.7,
                                       vacancy_type='V',
                                       inert_species=['Au', 'Pd'],
                                       verbose=True)
        self.assertEqual(mapped.get_chemical_formula(), 'H6Au4Pd4V2')


def suite():
    test_classes_to_run = [TestMapSites]
    suites_list = []
    for test_class in test_classes_to_run:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(test_class)
        suites_list.append(suite)
    test_suite = unittest.TestSuite(suites_list)
    return test_suite


if __name__ == '__main__':
    unittest.main()
