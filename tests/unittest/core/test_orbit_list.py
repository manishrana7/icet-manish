import unittest
from ase.build import bulk
from icet.core.orbit_list import create_orbit_list
from icet.core.cluster import Cluster
from icet import Structure

class TestOrbitList(unittest.TestCase):
    """
    Test class for testing orbit list
    """

    def __init__(self, *args, **kwargs):
        super(TestOrbitList, self).__init__(*args, **kwargs)
        self.atoms = bulk("Al")
        self.structure = Structure.from_atoms(self.atoms)
        self.cutoffs = [10,10]

    def setUp(self):
        """
        Setup test.
        """
        self.orbit_list = create_orbit_list(self.structure, self.cutoffs)

    def test_equivalent_sites_size(self):
        """
        Test that all the equivalent sites have the same geometrical size
        """

        for orbit in self.orbit_list.orbits:
            orbit.representative_cluster.geometrical_size
if __name__ == '__main__':
    unittest.main()
