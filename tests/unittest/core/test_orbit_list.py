import unittest
from ase.build import bulk
from icet.core.orbit_list import create_orbit_list
from icet.core.cluster import Cluster
from icet import Structure
from icet.tools.geometry import get_permutation
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
            size = orbit.geometrical_size
            for eq_sites in orbit.equivalent_sites:
                cluster = Cluster(self.structure, eq_sites,True)
                self.assertAlmostEqual(cluster.geometrical_size, size,places=5)
    
    def test_allowed_permutations(self):
        """
        Test allowed permutations of orbit.
        """

        for orbit in self.orbit_list.orbits:
            rep_sites = orbit.representative_sites
            translated_sites = self.orbit_list.get_sites_translated_to_unit_cell(rep_sites,False)
            permutations = orbit.allowed_permutations
            for perm in permutations:
                perm_sites = get_permutation(rep_sites, perm)
                self.assertIn(perm_sites, translated_sites)
                
if __name__ == '__main__':
    unittest.main()
