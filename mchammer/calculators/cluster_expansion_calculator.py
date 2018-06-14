from mchammer.calculators.base_calculator import BaseCalculator
from typing import List


class ClusterExpansionCalculator(BaseCalculator):
    """
    Cluster expansion calculator.

    Class for efficiently calculating the
    cluster expansion property
    for a specific structure

    Parameters
    ----------
    atoms : ASE Atoms object
        the structure that the calculator can use
        to optimize the calculate functions.

    cluster_expansion : icet ClusterExpansion object
    name : str
        human readable identifier for this calculator
    scaling : float (default len(atoms))
        this scales the property that is gotten from
        cluster vector times ECIs. By default the
        cluster vector times ECIs is assumed to give
        property/atom and thus the default value is
        multiplied by number of atoms.


    Todo
    ----
    * add the real occupation constraints when
      that is setup in the cluster space.

    """

    def __init__(self, atoms, cluster_expansion,
                 name='Cluster Expansion Calculator', scaling=None):
        super().__init__(atoms=atoms, name=name)
        self._cluster_expansion = cluster_expansion
        if scaling is None:
            self._property_scaling = len(atoms)
        else:
            self._property_scaling = scaling

    @property
    def cluster_expansion(self):
        """
        icet ClusterExpansion object.
        """
        return self._cluster_expansion

    def calculate_total(self, *, occupations: List[int]) -> float:
        """
        Calculates the total property of the current configuration.

        Parameters
        ----------
        occupations: list of int
            the entire occupation vector (i.e. list of atomic species)

        Returns
        -------
        total value of the property
        """
        self.atoms.set_atomic_numbers(occupations)
        return self.cluster_expansion.predict(self.atoms) * \
            self._property_scaling

    def calculate_local_contribution(self, local_indices: List[int] = None,
                                     occupations: List[int] = None) -> float:
        """
        Returns the sum of the contributions from the indices in the input
        list. `local_indices` refers to the lattice sites from which the local
        contributions should be summed up from. Occupations is the entire
        occupation vector.

        Parameters
        ----------
        local_indices : list of int
            the lattice indices for which to obtain the local contribution
        occupations : list of int
            the entire occupation vector

        Returns
        -------
        sum of contributions
        """
        if local_indices is None:
            raise TypeError("Missing required keyword argument: local_indices")
        if occupations is None:
            raise TypeError("Missing required keyword argument: occupations")
        return self.calculate_total(occupations=occupations) * \
            self._property_scaling

    @property
    def occupation_constraints(self):
        """A map from site to allowed species."""
        elements = list(
            self.cluster_expansion.cluster_space.element_map.keys())
        return [elements] * len(self.atoms)
