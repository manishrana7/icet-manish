import unittest
from mchammer.configuration_manager import ConfigurationManager
from ase.build import bulk
from ase.atoms import Atoms


class TestConfigurationManager(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestConfigurationManager, self).__init__(*args, **kwargs)
        self.atoms = bulk("Al").repeat(3)
        self.constraints

    def setUp(self):
        


