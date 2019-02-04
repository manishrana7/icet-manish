import unittest

from ase.build import bulk
from icet import ClusterSpace
from mchammer.observers import ClusterCountObserver


class TestClusterCountObserver(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestClusterCountObserver, self).__init__(*args, **kwargs)

        self.atoms = bulk('Al')

        cutoffs = [6, 6]
        subelements = ['Al', 'Ge']
        self.cs = ClusterSpace(self.atoms, cutoffs, subelements)
        self.tag = "count_dracula"
        self.interval = 10

    def shortDescription(self):
        """Silences unittest from printing the docstrings in test cases."""
        return None

    def setUp(self):
        """Set up observer before each test."""
        self.observer = ClusterCountObserver(
            cluster_space=self.cs, atoms=self.atoms, tag=self.tag, interval=self.interval)

    def test_property_tag(self):
        """Tests property tag."""
        self.assertEqual(self.observer.tag, self.tag)

    def test_property_interval(self):
        """Tests property interval."""
        self.assertEqual(self.observer.interval, self.interval)

    def test_get_observable(self):
        """Tests observable is returned accordingly."""
        counts = self.observer.get_observable(self.atoms)
        print("print observable")
        for tag, count in counts.items():
            print(tag, count)

        print("print observable")


if __name__ == '__main__':
    unittest.main()
