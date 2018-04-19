import unittest
from mchammer.configuration_manager import ConfigurationManager
from ase.build import bulk


class TestConfigurationManager(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestConfigurationManager, self).__init__(*args, **kwargs)
        self.atoms = bulk("Al").repeat(3)
        self.constraints = [[13, 47]*len(self.atoms)]
        self.strict_constraints = self.constraints
        self.sublattices = [[i for i in range(len(self.atoms))]]

    def setUp(self):
        self.cm = ConfigurationManager(
            self.atoms.numbers, self.strict_constraints, self.sublattices,
            self.constraints)

    def test_type(self):
        """ Test cm type """
        self.assertIsInstance(self.cm, ConfigurationManager)

    def test_check_occupation_constraint(self):
        """ Test the check occupation constraint method. """

        # Check that equal constraint is allowed
        constraint = [[1, 2], [1, 2], [1], [2]]
        strict_constraint = [[1, 2], [1, 2], [1], [2]]
        self.cm._check_occupation_constraint(strict_constraint, constraint)

        # Check that a more accepting constraint throws ValueError
        constraint = [[1, 2], [1, 2], [1], [2, 3]]
        strict_constraint = [[1, 2], [1, 2], [1], [2]]
        with self.assertRaises(Exception) as context:
            self.cm._check_occupation_constraint(strict_constraint, constraint)

        self.assertTrue("User defined occupation constraints must be "
                        "stricter or as strict as strict"
                        " occupations constraints." in str(context.exception))

        # Check that the length of the constraints throw a value error
        constraint = [[1, 2], [1, 2], [1]]
        strict_constraint = [[1, 2], [1, 2], [1], [2]]
        with self.assertRaises(Exception) as context:
            self.cm._check_occupation_constraint(strict_constraint, constraint)

        self.assertTrue("strict occupations and occupation "
                        "constraints must be equal length"
                        in str(context.exception))

    def test_property_occupations(self):
        """
        Test that the occupation property returns expected result
        and that it cannot be modified.
        """

        self.assertListEqual(list(self.cm.occupations),
                             list(self.atoms.numbers))

        # Test that the property can not be set.
        with self.assertRaises(AttributeError) as context:
            self.cm.occupations = []
        self.assertTrue("can't set attribute" in str(context.exception))

        # Test that the property can't be set by modifying the
        # returned occupations
        occupations = self.cm.occupations
        occupations[0] = occupations[0] + 1
        self.assertNotEqual(list(occupations), list(self.cm.occupations))


if __name__ == '__main__':
    unittest.main()
