"""
This test checks that a cluster expansion model can be initialized
with any structure in the test database and can predict a property.
"""
import unittest
import tempfile

from ase.db import connect
from icet import ClusterSpace, ClusterExpansion
from ase.build import bulk
from io import StringIO


def strip_surrounding_spaces(input_string):
    """
    Helper function that removes both leading and trailing spaces from a
    multi-line string.

    Returns
    -------
    str
        original string minus surrounding spaces and empty lines
    """
    s = []
    for line in StringIO(input_string):
        if len(line.strip()) == 0:
            continue
        s += [line.strip()]
    return '\n'.join(s)


class TestClusterExpansion(unittest.TestCase):
    """
    Container for tests of the class functionality
    """

    def __init__(self, *args, **kwargs):
        super(TestClusterExpansion, self).__init__(*args, **kwargs)
        self.db = connect('structures_for_testing.db')
        self.subelements = ['H', 'He', 'Pb']
        self.cutoffs = [1.4] * 3

    def test_init(self):
        """Test that initialization works."""
        atoms = bulk('Au')
        cutoffs = [3.0]
        cs = ClusterSpace(atoms, cutoffs, ['Au', 'Pd'])
        with self.assertRaises(ValueError) as context:
            ClusterExpansion(cs, [0.0])
        msg = 'cluster_space and parameters must have the same' + \
            ' length (3 != 1)'
        self.assertEqual(str(context.exception), msg)

    def test_all_systems_in_database(self):
        """
        Test all structures in database
        """
        for row in self.db.select():
            atoms_row = row.toatoms()
            atoms_row.set_chemical_symbols(
                len(atoms_row) * [self.subelements[0]])
            self._test_clusterexpansion_model(
                atoms_row, self.cutoffs, self.subelements)

    def _test_clusterexpansion_model(self, atoms, cutoffs, subelements):
        """
        Test cluster expansion init and prediction

        Parameters
        ----------
        atoms : ASE Atoms object
            atomic configuration
        subelements : list of strings (chemical symbols)
            list of elements that are allowed
        """
        cs = ClusterSpace(atoms, cutoffs, subelements)
        params_len = cs.get_cluster_space_size()
        params = list(range(params_len))

        ce = ClusterExpansion(cs, params)
        predicted_val = ce.predict(atoms)

        self.assertIsInstance(predicted_val, float)

    def test_read_write(self):
        """
        Test read write functionality.
        """
        atoms = bulk("Al")
        cluster_space = ClusterSpace(atoms, self.cutoffs, self.subelements)
        params = list(range(len(cluster_space)))
        ce = ClusterExpansion(cluster_space, params)

        f = tempfile.NamedTemporaryFile()
        ce.write(f.name)
        f.seek(0)
        ce_read = ClusterExpansion.read(f.name)

        # Test clusterspace in cluster expansion
        self.assertEqual(cluster_space._atoms,
                         ce_read.cluster_space._atoms)
        self.assertEqual(list(cluster_space._cutoffs),
                         list(ce_read.cluster_space._cutoffs))
        self.assertEqual(cluster_space._chemical_symbols,
                         ce_read.cluster_space._chemical_symbols)
        self.assertEqual(cluster_space._mi, ce_read.cluster_space._mi)
        self.assertEqual(cluster_space._verbosity,
                         ce_read.cluster_space._verbosity)

        # Test equal parameters
        self.assertEqual(ce_read.parameters, params)

    def test_repr(self):
        """
        Testing repr functionality
        """
        atoms = bulk("Al")
        cutoffs = [3.0] * 3
        subelements = ['Al', 'Pd']
        cluster_space = ClusterSpace(atoms, cutoffs, subelements)
        params = list(range(len(cluster_space)))
        ce = ClusterExpansion(cluster_space, params)

        retval = ce.__repr__()
        target = """
================================= Cluster Expansion =================================
 chemical species: Al Pd
 cutoffs: 3.0000 3.0000 3.0000
 total number of orbits: 5
 number of orbits by order: 0= 1  1= 1  2= 1  3= 1  4= 1
-------------------------------------------------------------------------------------
index | order |  radius  | multiplicity | orbit_index | multi_component_vector | ECI 
-------------------------------------------------------------------------------------
   0  |   0   |   0.0000 |        1     |      -1     |           .            |  0  
   1  |   1   |   0.0000 |        1     |       0     |          [0]           |  1  
   2  |   2   |   1.4319 |        6     |       1     |         [0, 0]         |  2  
   3  |   3   |   1.6534 |        8     |       2     |       [0, 0, 0]        |  3  
   4  |   4   |   1.7537 |        2     |       3     |      [0, 0, 0, 0]      |  4  
=====================================================================================
""" # noqa
        self.assertEqual(strip_surrounding_spaces(target),
                         strip_surrounding_spaces(retval))

    def test_get_string_representation(self):
        """
        Testing _get_string_representation functionality.
        """
        atoms = bulk("Al")
        cutoffs = [3.0] * 3
        subelements = ['Al', 'Pd']
        cluster_space = ClusterSpace(atoms, cutoffs, subelements)
        params = list(range(len(cluster_space)))
        ce = ClusterExpansion(cluster_space, params)

        retval = ce._get_string_representation(print_threshold=2,
                                               print_minimum=1)
        target = """
================================= Cluster Expansion =================================
 chemical species: Al Pd
 cutoffs: 3.0000 3.0000 3.0000
 total number of orbits: 5
 number of orbits by order: 0= 1  1= 1  2= 1  3= 1  4= 1
-------------------------------------------------------------------------------------
index | order |  radius  | multiplicity | orbit_index | multi_component_vector | ECI 
-------------------------------------------------------------------------------------
   0  |   0   |   0.0000 |        1     |      -1     |           .            |  0  
 ...
   4  |   4   |   1.7537 |        2     |       3     |      [0, 0, 0, 0]      |  4  
=====================================================================================
""" # noqa
        self.assertEqual(strip_surrounding_spaces(target),
                         strip_surrounding_spaces(retval))


if __name__ == '__main__':
    unittest.main()
