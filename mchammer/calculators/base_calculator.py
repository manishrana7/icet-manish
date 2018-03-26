from abc import ABC, abstractmethod


class BaseCalculator(ABC):
    """
    Base ensemble abstract class.

    Attributes
    ----------
    name : str
        name of the calculator.    
    """

    def __init__(self, atoms, name='BaseCalculator'):
        self._atoms = atoms
        self.name = name
    
    @property
    def atoms(self):
        """
        ASE Atoms object
        """
        return self._atoms

    @abstractmethod
    def calculate_total(self):
        pass

    @abstractmethod
    def calculate_local_contribution(self):
        pass

    def set_elements(self, indices, elements):
        """
        Set elements on the atoms object.

        indices : list of int
            The list refer to indices on the lattice.
        elements : list of int
            The elements to put on the corresponding indice.
        """

        assert len(indices) == len(
            elements), "indices and elements need to be the same size."

        for index, element in zip(indices, elements):
            self.atoms.numbers[index] = element
