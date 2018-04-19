import unittest
from mchammer.configuration_manager import ConfigurationManager
from ase.build import bulk
from ase.atoms import Atoms


class TestConfigurationManager(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestConfigurationManager, self).__init__(*args, **kwargs)
        self.atoms = bulk("Al").repeat(3)
        self.constraints = [[13,47 ]*len(self.atoms)]
        self.strict_constraints = self.constraints
        self.sublattices = [[i for i in range(len(self.atoms))]]

    def setUp(self):
        self.cm = ConfigurationManager(self.atoms.numbers, self.strict_constraints, self.sublattices, self.constraints)

    def test_type(self):
        """ Test cm type """
        self.assertIsInstance(self.cm, ConfigurationManager)


if __name__ == '__main__':
    unittest.main()
