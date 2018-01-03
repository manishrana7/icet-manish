#!/usr/bin/env python3

import unittest

from icetdev.tools import ConvexHull
import numpy as np

concentrations_1 = [0.0, 1.0, 0.4, 0.6]
energies_1 = [0.0, 10.0, -10.0, 20.0]

concentrations_2 = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.1, 0.1], [0.3, 0.3]]
energies_2 = [0.0, 10.0, -10.0, 3.0, -7.0]


class TestConvexHull(unittest.TestCase):
    '''
    Container for tests of the class functionality
    '''

    def setUp(self):
        '''
        Instantiate class before each test.
        '''
        self.ch = ConvexHull(concentrations_1, energies_1)

    def test_init(self):
        '''
        Testing that the setup
        (initialization) of tested class work
        '''
        ch = ConvexHull(concentrations_1, energies_1)
        self.assertEqual(ch.dimensions, 1)
        self.assertTrue(np.allclose(ch.energies,
                                    np.array([0.0, -10.0, 10.0])))

        ch = ConvexHull(concentrations_2, energies_2)
        self.assertEqual(ch.dimensions, 2)
        self.assertTrue(np.allclose(ch.energies,
                                    np.array([0.0, 10.0, -10.0, -7.0])))

    def test_get_energy_at_convex_hull(self):
        '''
        Testing energy at convex hull retrieval functionality.
        '''
        convex_hull_energies = self.ch.get_energy_at_convex_hull([0.2, 0.7])
        self.assertTrue(np.allclose(convex_hull_energies,
                                    np.array([-5.0, 0.0])))


def suite():
    suites_list = []
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestConvexHull)
    suites_list.append(suite)
    test_suite = unittest.TestSuite(suites_list)
    return test_suite


if __name__ == '__main__':
    unittest.main()
