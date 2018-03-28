#!/usr/bin/env python3

import unittest

from mchammer import BaseEnsemble, DataContainer
from icet import ClusterSpace, ClusterExpansion
from ase.build import bulk
from collections import OrderedDict
from datetime import date


class TestDataContainer(unittest.TestCase):
    """
    Container for the tests of the class functionality

    """
    def __init__(self, *args, **kwargs):
        super(TestDataContainer, self).__init__(*args, **kwargs)
        subelements = ['W', 'I', 'Li']
        cutoffs = [1.4] * 3
        atoms_prim = bulk("Al")
        cluster_space = ClusterSpace(atoms_prim, cutoffs, subelements)
        params = list(range(len(cluster_space)))

        self.observers = ['some_observer']
        self.obs_tag = 'energy'
        self.obs_observable = -12.0

        self.ce = ClusterExpansion(cluster_space, params)
        # self.atoms = self.ce.primitive_structure.repeat(4)
        self.atoms = atoms_prim.repeat(4)

    def setUp(self):
        """
        setUp all tests with this
        """
        SomeEnsemble = BaseEnsemble(atoms=self.atoms,
                                    calculator=self.ce)
        self.data = DataContainer(SomeEnsemble)

    def test_add_structure(self):
        """
        Test that reference structure is added to DataContainer.
        """
        self.data.add_structure(self.atoms)
        self.assertEqual(self.data.structure, self.atoms)
        # test another type fails
        with self.assertRaises(Exception):
            self.data.add_structure('something')

    def test_add_observable(self):
        """
        Test that observables are added to DataContainer.
        """
        self.data.add_observable('energy', float)
        self.data.add_observable('temperature', float)
        self.assertEqual(len(self.data.observables), 2)

    def test_add_parameter(self):
        """
        Test that parameters are added to DataContainer
        """
        self.data.add_parameter('temperature', 1000.0)
        self.data.add_parameter('chemical-potential-difference', -0.5)
        self.assertEqual(len(self.data.parameters), 2)

    def test_add_metadata(self):
        """
        Tests that metadata is added to DataContainer.
        """
        today = date.today()
        author, hostname, time = 'user', 'host', today

        self.data.add_metadata('author', author)
        self.data.add_metadata('hostname', hostname)
        self.data.add_metadata('date', time)

        self.assertEqual(len(self.data.metadata), 3)

    def test_append_data(self):
        """
        Test that data is appended to DataContainer.
        """
        minimum_interval = 10
        obs_interval = 100
        for mcstep in range(1, 1001):
            if mcstep % minimum_interval == 0:
                for obs in self.observers:
                    if mcstep % obs_interval == 0:
                        self.data.append(mcstep,
                                         self.obs_tag,
                                         self.obs_observable)
                        self.data.update()

        self.assertEqual(len(self.data), minimum_interval)

    def test_parameters(self):
        """
        Test that added parameters has OrderedDict type.
        """
        target = OrderedDict([('temperature', 300.0),
                              ('chemical-potential-difference', -0.5)])

        self.data.add_parameter('temperature', 300.0)
        self.data.add_parameter('chemical-potential-difference', -0.5)

        retval = self.data.parameters
        self.assertEqual(retval, target)

    def test_observables(self):
        """
        Test that added observables has OrderedDict type.

        Todo
        ----
        * Add seed or PRNG, interval, minimum interval
        """
        target = OrderedDict([('energy', float),
                              ('temperature', float),
                              ('sro', list)])

        self.data.add_observable('energy', float)
        self.data.add_observable('temperature', float)
        self.data.add_observable('sro', list)
        retval = self.data.observables

        self.assertEqual(retval, target)

    def test_metadata(self):
        """
        Test added metadata.
        """
        target = OrderedDict([('author', 'user'),
                              ('hostname', 'host'),
                              ('date', date(2018, 3, 23))])

        this_date = date(2018, 3, 23)
        author, hostname, time = 'user', 'host', this_date

        self.data.add_metadata('author', author)
        self.data.add_metadata('hostname', hostname)
        self.data.add_metadata('date', time)
        retval = self.data.metadata

        self.assertDictEqual(retval, target)

    def test_get_data(self):
        """
        Test that data is returned from DataContainer
        has right format (list of lists).
        """
        target = [[100.0, -12.0], [200.0, -12.0], [300.0, -12.0],
                  [400.0, -12.0], [500.0, -12.0], [600.0, -12.0],
                  [700.0, -12.0], [800.0, -12.0], [900.0, -12.0]]

        minimum_interval = 10
        obs_interval = 100
        for mcstep in range(1, 1000):
            if mcstep % minimum_interval == 0:
                for obs in self.observers:
                    if mcstep % obs_interval == 0:
                        self.data.append(mcstep,
                                         self.obs_tag,
                                         self.obs_observable)
                        self.data.update()
        retval = self.data.get_data(['mcstep', self.obs_tag])
        self.assertListEqual(target, retval)

        with self.assertRaises(AssertionError):
            self.data.get_data(['temperature'])

    def test_reset(self):
        """
        Test this functionality removes the
        appended data to DataContainer.
        """
        temperatures = list(range(1000, 100, -100))
        for mcstep, temperature in enumerate(temperatures):
            self.data.append(mcstep, 'temperature', temperature)
        self.data.reset()
        self.assertEqual(len(self.data), 0)

    def test_len(self):
        """
        Test total number of rows is returned.
        """
        temperatures = list(range(1000, 100, -100))
        for mcstep, temperature in enumerate(temperatures):
            self.data.append(mcstep, 'temperature', temperature)
            self.data.update()
        self.assertEqual(len(self.data), len(temperatures))

    def test_get_average(self):
        """
        Test functionality.
        """
        pass

    def test_load_data(self):
        """
        Test read data container from file.
        """
        pass

    def test_save_data(self):
        """
        Test data container is saved in file.
        """
        pass


if __name__ == '__main__':
    unittest.main()
