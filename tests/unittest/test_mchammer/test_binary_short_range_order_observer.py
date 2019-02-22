import unittest

from ase.build import bulk
from icet import ClusterSpace
from mchammer.observers import BinaryShortRangeOrderObserver


class TestBinaryShortRangeOrderObserver(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestBinaryShortRangeOrderObserver,
              self).__init__(*args, **kwargs)

        self.atoms = bulk('Au').repeat([2, 1, 1])
        self.atoms[1].symbol = 'H'
        self.atoms = self.atoms.repeat(2)

        cutoffs = [3]
        subelements = [['Au', 'Pd'], ['H']]*(len(self.atoms)//2)
        self.cs = ClusterSpace(self.atoms, cutoffs, subelements)
        self.interval = 10
        self.radius = 5

    def shortDescription(self):
        """Silences unittest from printing the docstrings in test cases."""
        return None

    def setUp(self):
        """Set up observer before each test."""
        self.observer = BinaryShortRangeOrderObserver(
            cluster_space=self.cs, structure=self.atoms,
            interval=self.interval, radius=self.radius)

    def test_property_tag(self):
        """Tests property tag."""
        self.assertEqual(self.observer.tag, "BinaryShortRangeOrderObserver")

    def test_property_interval(self):
        """Tests property interval."""
        self.assertEqual(self.observer.interval, self.interval)

    def test_get_observable(self):
        """Tests observable is returned accordingly."""
        sro_parameters = self.observer.get_observable(self.atoms)

        self.assertEqual(sro_parameters['sro_Au_1'], 1.0)
        self.assertEqual(sro_parameters['sro_Au_2'], 1.0)
        self.assertEqual(sro_parameters['sro_Au_3'], 1.0)

    def test_get_concentrations(self):
        conc = self.observer._get_concentrations(self.atoms)

        self.assertEqual(conc['Au'], 1.0)
        self.assertEqual(conc['Pd'], 0.0)
        self.assertEqual(conc['H'], 1.0)


if __name__ == '__main__':
    unittest.main()
