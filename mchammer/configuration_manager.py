import copy
import random
from ase import Atoms
from typing import Dict, List, Tuple, Union


class SwapNotPossibleError(Exception):
    pass


class ConfigurationManager(object):
    """The ConfigurationManager store its own state of the
    configuration that is being sampled in the ensemble.

    ConfigurationManager is responsible for all information and
    handling of the configuration.

    Parameters
    ----------
    atoms : ASE Atoms
        configuration to be handled
    strict_constraints : list of list of int
        strictest form of the allowed occupations
    sublattices : list of list of int
        the integers refer to lattice site indices in the configuration.
    occupation_constraints : list of list of int
        optional occupation constraint to enfore a more stricter species
        occupation than what is allowed from the Calculator.

    Todo
    ----
    * occupation constraint not implemented
    * add check that all sites in the different sublattices all have the same
      occupation constraint.
    * swap "element" for "species"
    * revise docstrings
    * clarify "sublattices" and "occupation_constraints";
      the OccupationConstraints class should help here
    """

    def __init__(self, atoms: Atoms,
                 strict_constraints: Union[List[list], List[int]],
                 sublattices: Union[List[list], List[int]],
                 occupation_constraints: List[List[int]]=None):

        self._atoms = atoms
        self._occupations = atoms.numbers
        self._sublattices = sublattices

        if occupation_constraints is not None:
            self._check_occupation_constraint(
                strict_constraints, occupation_constraints)
        else:
            occupation_constraints = strict_constraints
        self._occupation_constraints = occupation_constraints
        self._possible_elements = self._set_up_possible_elements()
        self._element_occupation = self._get_element_occupation()

    def _set_up_possible_elements(self)->List[int]:
        """
        Returns a list of the possible elements in
        the entire configuration.
        """
        possible_element = set()
        for occ in self._occupation_constraints:
            for element in occ:
                possible_element.add(element)

        return list(possible_element)

    def _get_element_occupation(self) -> List[Dict[int, int]]:
        """
        Return the element occupation, i.e. the element -> lattice sites dict.
        """
        element_occupations = []
        for i, sublattice in enumerate(self._sublattices):
            element_dict = {}  # type: dict
            for element in self._possible_elements:
                element_dict[element] = []
            for lattice_site in sublattice:
                element = self._occupations[lattice_site]
                element_dict[element].append(lattice_site)
            element_occupations.append(element_dict)
        return element_occupations

    def _check_occupation_constraint(self, strict_constraints,
                                     occupation_constraints):
        """
        Checks that the user defined occupations constrains are stricter or
        as strict as the strict occupations.

        Parameters
        ----------
        strict_constraints : list of list of int
        occupation_constraints : list of list of int
        """

        if not len(strict_constraints) == len(occupation_constraints):
            raise ValueError(
                "strict occupations and occupation"
                " constraints must be equal length")

        for strict_occ, occ in zip(strict_constraints, occupation_constraints):
            if not set(occ).issubset(strict_occ):
                raise ValueError(
                    "User defined occupation constraints must be "
                    "stricter or as strict as strict occupations constraints.")

    @property
    def occupations(self) ->List[int]:
        """The occupations of the configuration."""
        return self._occupations.copy()

    @property
    def occupation_constraints(self) ->List[List[int]]:
        """The occupation constraints for this configuration manager."""
        return copy.deepcopy(self._occupation_constraints)

    @property
    def sublattices(self) -> List[List[int]]:
        """The sublattices of the configuration."""
        return copy.deepcopy(self._sublattices)

    @property
    def atoms(self):
        """The atoms object."""
        atoms = self._atoms.copy()
        atoms.set_atomic_numbers(self.occupations)
        return atoms

    def get_swapped_state(self, sublattice: int) -> Tuple[List[int],
                                                          List[int]]:
        """
        Returns a tuple of a list of two random indices in a specific
        sublattice and what elements they will occupy after a swap.
        The two indices refer to lattice sites and the indices
        will, if swapped, produce a new, different configuration
        that is allowed by the occupation constraints.

        Parameters
        ----------
        sublattice : int
            The sublattice to pick indices from.
        """
        try:
            index_1 = random.choice(self._sublattices[sublattice])
        except IndexError:
            raise SwapNotPossibleError("Sublattice is empty")
        # assuming all sites in this sublattice have same allowed occupations

        possible_swap_elements = set(
            self._occupation_constraints[index_1]) - set(
                [self._occupations[index_1]])

        total_number_of_indices = 0
        for element in possible_swap_elements:
            total_number_of_indices += len(
                self._element_occupation[sublattice][element])
        try:
            index = random.randint(0, total_number_of_indices - 1)
        except ValueError:
            raise SwapNotPossibleError
        for element in possible_swap_elements:
            if index < len(self._element_occupation[sublattice][element]):
                index_2 = self._element_occupation[sublattice][element][index]
                break
            else:
                index -= len(self._element_occupation[sublattice][element])
        else:
            raise SwapNotPossibleError
        return [index_1, index_2], [self._occupations[index_2],
                                    self._occupations[index_1]]

    def get_flip_state(self, sublattice: int) -> Tuple[int, int]:
        """
        Returns a tuple of a list of a random index in a specific sublattice
        and a list of the element to flip it to.

        Parameters
        ----------
        sublattice
            sublattice from which to pick indices

        Todo
        ----
        * improve description of "list of the element to flip it to"
        """

        index = random.choice(self._sublattices[sublattice])
        element = random.choice(list(
            set(self._occupation_constraints[index]) - set(
                [self._occupations[index]])))
        return index, element

    def update_occupations(self, list_of_sites: List[int],
                           list_of_elements: List[int]):
        """
        Updates the element occupation of the configuration being sampled.
        This will change the state in both the configuration in the calculator
        and the configuration manager.

        Parameters
        ----------
        list_of_sites
            list of indices to be changed in the configuration
        list_of_elements
            list of species to be assigned to lattice sites

        Todo
        ----
        * change occupation
        * change the element occupation
        """

        # change element occupation dict
        for index, element in zip(list_of_sites, list_of_elements):
            current_element = self._occupations[index]
            for sublattice_index, sublattice in enumerate(self._sublattices):
                if index in sublattice:
                    # Remove index from current element
                    self._element_occupation[
                        sublattice_index][current_element].remove(
                        index)
                    # Add index to new element
                    self._element_occupation[sublattice_index][element].append(
                        index)

        # change occupation list
        for index, element in zip(list_of_sites, list_of_elements):
            self._occupations[index] = element
