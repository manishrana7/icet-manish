#!/usr/bin/env python3

import unittest

from mchammer import BaseEnsemble, DataContainer
from icet import ClusterSpace, ClusterExpansion
from ase.build import bulk


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

        self.ce = ClusterExpansion(cluster_space, params)
        # self.atoms = self.ce.primitive_structure.repeat(4)
        self.atoms = atoms_prim.repeat(4)

    def setUp(self):
        """
        setUp all tests with this
        """
        SomeEnsemble = BaseEnsemble(atoms=self.atoms,
                                    calculator=self.ce,
                                    temperature=100.0,
                                    interval=20)
        self.data_container = DataContainer(SomeEnsemble)

    def test_add_structure(self):
        """
        Test functionality.
        """
        pass

    def test_add_observable(self):
        """
        Test functionality.
        """
        pass

    def test_add_parameter(self):
        """
        Test functionality.
        """
        pass

    def test_add_metadata(self):
        """
        Test that metadata (date, machine, ensembles)
        is stored here properly.
        """
        pass

    def test_append_data(self):
        """
        Test functionality.
        """
        pass

    def test_parameters(self):
        """
        Test functionality.
        """
        pass

    def test_observers(self):
        """
        Test functionality.
        """
        pass

    def test_get_data(self):
        """
        Test functionality.
        """
        pass

    def test_get_average(self):
        """
        Test functionality.
        """
        pass

    def test_reset_data_container(self):
        """
        Test functionlity.
        """
        pass

    def test_load_data_container(self):
        """
        Test read data container from file.
        """
        pass

    def test_save_data_container(self):
        """
        Test data container is saved in file.
        """
        pass


if __name__ == '__main__':
    unittest.main()
