import random


class ConfigurationManager(object):
    """
    ConfigurationManager store its own state of the
    configuration that are being sampled in the ensemble.

    ConfigurationManager is responsible for all information and
    handling of the configuration.

    Parameters
    ----------
    occupations : list of int
        The occupation of the configurations. The integers
        refer to elements via their number in the periodic table
    strict_occupations : list of list of int
        the strictest form of the allowed occupations.
    occupation_constraints : list of list of int
        optional occupation constraint to enfore a more stricter species
        occupation than what is allowed from the Calculator.
    sublattices : list of list of int
        the integers refer to lattice site indices in the configuration.

    Todo
    ----
    * occupation constraint not implemented
    * add check that all sites in the different sublattices all have the same
      occupation constraint.
    """

        def __init__(self, occupations, strict_occupations, sublattices,  occupation_constraints=None):

            self._occupations = occupations
            self._sublattices = sublattices

            if occupation_constraints is not None:
                self._check_occupation_constraint(
                    strict_occupations, occupation_constraints)
            self._occupation_constraints = occupation_constraints
            self._possible_elements = self._get_possible_elements()
            self._element_occupation = self._get_element_occupation()

            # internal representation
            # (strings are probably not preferable in actual implementation;
            # they are used here for clarity)
            # ConfigurationManager._sublattices = [['Pd', 'Au'],
            #                                      ['H',  'X']]
            # ConfigurationManager._sites = {'Pd': [0, 1, 2, 3],
            #                                'Au': [4, 5, 6, 7],
            #                                'H':  [8, 9, 10],
            #                                'X':  [11, 12, 13, 14]}

        def _get_possible_elements(self):
            """
            Return a list of the possible elements in
            the entire configuration.
            """
            possible_element = set()
            for occ in self.occupation_constraints:
                for element in occ:
                    possible_element.insert(element)

            return list(possible_element)

        def _get_element_occupation(self):
            """
            Return the element occupation, i.e. the element -> lattice sites dict.
            """
            element_occupations = []
            for i, sublattice in enumerate(self.sublattices)
                element_dict = {}
                for element in self._possible_elements:
                    element_dict[element] = []
                for lattice_site in sublattice:
                    element = self.occupations[lattice_site]
                    element_dict[element].push_back(lattice_site)
                element_occupations.push_back(element_dict)
            return element_occupations
        def _check_occupation_constraint(self, strict_occupations, occupation_constraints):
            """
            Checks that the user defined occupations constrains are stricter or
            as strict as the strict occupations.

            Parameters
            ----------
            strict_occupations : list of list of int
            occupation_constraints : list of list of int
            """
            error_message = "strict occupations and occupation constraints must be equal length"
            assert len(strict_occupations) == len(
                occupation_constraints), error_message

            for strict_occ, occ in zip(strict_occupations, occupation_constraints):
                if not set(occ).issubset(strict_occ):
                    raise RuntTimeEerror(
                        "User defined occupation constraints must be "
                        "stricter or as strict as strict occupations constraints.")

        @property
        def occupations(self):
            """ The occupations of the configuration."""
            return self._occupations.copy()

        @property
        def occupation_constraints(self):
            """ The occupation constraints for this configuration manager."""
            return self._occupation_constraints.copy()

        @property
        def sublattices(self):
            """ The sublattices of the configuration. """
            return self._sublattices.copy()

        def get_swap_indices(self, sublattice: int):
            """
            Returns a tuple of two random indices in a specific sublattice.

            Parameters
            ----------
            sublattice : int
                The sublattice to pick indices from.
            """

            index_1 = np.random.choice(self.sublattices[sublattice])
            element_1 = self.occupations[index_1]
            # assuming all sites in this sublattice have same allowed occupations

            possible_swap_elements = set(self.occupation_constraints[index_1]) - set(self.occupations[index_1]))
            total_number_of_indices = 0
            for element in possible_swap_elements:
                total_number_of_indices += len(self._element_occupation[sublattice][element])
            index = random.randint(0, total_number_of_indices-1)
            for element in possible_swap_elements:
                if index < len(self._element_occupation[sublattice][element]):
                    index_2 = self._element_occupation[sublattice][element][index]
                    break
                else:
                    index -= len(self._element_occupation[sublattice][element])
            else:
                raise SwapNotPossibleError
            return index_1, index_2
            
                
            self._element_occupation[sublattice][]

            raise NotImplementedError

        def get_flip_index(self, sublattice: int):
            """
            Returns a tuple of a random index in a specific sublattice
            and the element to flip it to.

            Parameters
            ----------
            sublattice : int
                The sublattice to pick indices from.

            return : (int, int) (index, element)

            """

            index=np.random.choice(self.sublattices[sublattice])
            element=np.random.choice(
                set(self.occupation_constraints[index]) - set(self.occupations[index]))
            return index, element
