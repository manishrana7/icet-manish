import unittest

from ase.build import bulk
from icet import ClusterSpace
from mchammer.observers import ClusterCountObserver


class TestClusterCountObserver(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestClusterCountObserver, self).__init__(*args, **kwargs)

        self.atoms = bulk('Al').repeat([2, 1, 1])
        self.atoms[1].symbol = 'Ge'

        cutoffs = [3]
        subelements = ['Al', 'Ge', 'Si']
        self.cs = ClusterSpace(self.atoms, cutoffs, subelements)
        self.interval = 10

    def shortDescription(self):
        """Silences unittest from printing the docstrings in test cases."""
        return None

    def setUp(self):
        """Set up observer before each test."""
        self.observer = ClusterCountObserver(
            cluster_space=self.cs, atoms=self.atoms, interval=self.interval)

    def test_property_tag(self):
        """Tests property tag."""
        self.assertEqual(self.observer.tag, "ClusterCountObserver")

    def test_property_interval(self):
        """Tests property interval."""
        self.assertEqual(self.observer.interval, self.interval)

    def test_get_observable(self):
        """Tests observable is returned accordingly."""
        self.observer.get_observable(self.atoms)


if __name__ == '__main__':
    unittest.main()
