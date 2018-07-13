from ase import Atoms
from icet import ClusterExpansion
from mchammer.calculators.base_calculator import BaseCalculator
from typing import List


class ClusterExpansionCalculator(BaseCalculator):
    """A ClusterExpansionCalculator object enables the efficient
    calculation of properties described by a cluster expansion. It is
    specific for a particular (supercell) structure and commonly
    employed when setting up a Monte Carlo simulation, see
    :ref:`ensembles`.

    Cluster expansions, e.g., of the energy, typically yield property
    values *per site*. When running a Monte Carlo simulation one,
    however, considers changes in the *total* energy of the
    system. The default behavior is therefore to multiply the output
    of the cluster expansion by the number of sites. This behavior can
    be changed via the ``scaling`` keyword parameter.

    Parameters
    ----------
    atoms : :class:`ase:Atoms`
        structure for which to set up the calculator
    cluster_expansion : ClusterExpansion
        cluster expansion from which to build calculator
    name : str
        human-readable identifier for this calculator
    scaling : float
        scaling factor applied to the property value predicted by the
        cluster expansion

    Todo
    ----
    * add OccupationConstraints once available

    """

    def __init__(self, atoms: Atoms, cluster_expansion: ClusterExpansion,
                 name: str='Cluster Expansion Calculator',
                 scaling: float=None):
        super().__init__(atoms=atoms, name=name)
        self._cluster_expansion = cluster_expansion
        if scaling is None:
            self._property_scaling = len(atoms)
        else:
            self._property_scaling = scaling

    @property
    def cluster_expansion(self) -> ClusterExpansion:
        """ cluster expansion from which calculator was constructed """
        return self._cluster_expansion

    def calculate_total(self, *, occupations: List[int]) -> float:
        """
        Calculates and returns the total property value of the current
        configuration.

        Parameters
        ----------
        occupations
            the entire occupation vector (i.e. list of atomic species)
        """
        self.atoms.set_atomic_numbers(occupations)
        return self.cluster_expansion.predict(self.atoms) * \
            self._property_scaling

    def calculate_local_contribution(self, local_indices: List[int] = None,
                                     occupations: List[int] = None) -> float:
        """
        Calculates and returns the sum of the contributions to the property
        due to the sites specified in `local_indices`

        Parameters
        ----------
        local_indices
            sites over which to sum up the local contribution
        occupations
            entire occupation vector
        """
        if local_indices is None:
            raise TypeError('Missing required argument: local_indices')
        if occupations is None:
            raise TypeError('Missing required argument: occupations')
        return self.calculate_total(occupations=occupations) * \
            self._property_scaling

    @property
    def occupation_constraints(self) -> List[List[int]]:
        """ map from site to allowed species """
        species = list(
            self.cluster_expansion.cluster_space.species_map.keys())
        return [species] * len(self.atoms)
