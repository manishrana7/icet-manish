"""
Definition of the variance-constrained semi-grand canonical ensemble class.
"""

import numpy as np

from ase import Atoms
from ase.data import atomic_numbers, chemical_symbols
from ase.units import kB
from typing import List, Dict

from .. import DataContainer
from .base_ensemble import BaseEnsemble
from ..calculators.base_calculator import BaseCalculator


class VCSGCEnsemble(BaseEnsemble):
    """Variance-constrained semi-grand canonical (VCSGC) ensemble.

    Instances of this class allow one to simulate systems in the VCSGC
    ensemble (:math:`N\phi\kappa VT`), i.e. at constant temperature
    (:math:`T`), total number of sites (:math:`N=\sum_i N_i`),
    and two additional parameters :math:`\phi` and :math:`\kappa`, which
    constrain the concentration and variance of the concentration,
    respectively. The derivative of the canonical free energy can be
    expressed in observables of the ensemble,

    .. math::

        \\frac{1}{N} \\frac{\\partial F}{\\partial c} = - \phi - 2 N \kappa
        \\langle c \\rangle.

    Unlike the SGC ensemble, the VCSGC ensemble allows for sampling across
    multi-phase regions, meaning that the free energy can be recovered by
    direct integration even when such regions are present.

    The VCSGC ensemble currently supports systems with no more than two
    different species.

    When using this ensemble, please cite 
    Sadigh, B. and Erhart, P., Phys. Rev. B **86**, 134204 (2012)
    [SadErh12]_.

    Attributes
    -----------
    temperature : float
        temperature :math:`T` in appropriate units [commonly Kelvin]
    boltzmann_constant : float
        Boltzmann constant :math:`k_B` in appropriate
        units, i.e. units that are consistent
        with the underlying cluster expansion
        and the temperature units [default: eV/K]
    variance_parameter : float
        parameter that constrains the fluctuations of the concentration
    """

    def __init__(self, atoms: Atoms=None, calculator: BaseCalculator=None,
                 name: str='Variance-constrained semi-grand' +
                           ' canonical ensemble',
                 data_container: DataContainer=None, random_seed: int=None,
                 data_container_write_period: float=np.inf,
                 ensemble_data_write_interval: int=None,
                 trajectory_write_interval: int=None,
                 boltzmann_constant: float=kB, *, temperature: float,
                 concentration_parameters: Dict[str, float],
                 variance_parameter: float):

        super().__init__(
            atoms=atoms, calculator=calculator, name=name,
            data_container=data_container,
            random_seed=random_seed,
            data_container_write_period=data_container_write_period,
            ensemble_data_write_interval=ensemble_data_write_interval,
            trajectory_write_interval=trajectory_write_interval)

        self.temperature = temperature
        self.boltzmann_constant = boltzmann_constant

        self._concentration_parameters = None
        self.concentration_parameters = concentration_parameters
        self.variance_parameter = variance_parameter

        if len(self.configuration._allowed_species) > 2:
            raise NotImplementedError('VCSGCEnsemble does not yet support '
                                      'cluster spaces with more than two '
                                      'species.')
        self._species_counts = None
        self.species_counts = self.configuration.occupations

    def _do_trial_step(self):
        """ Carries out one Monte Carlo trial step. """
        self.total_trials += 1

        # choose swap
        sublattice_index = self.get_random_sublattice_index()
        index, new_species = \
            self.configuration.get_flip_state(sublattice_index)
        old_species = self.configuration.occupations[index]

        # Calculate difference in VCSGC thermodynamic potential.
        # Note that this assumes that only one atom was flipped.
        N = len(self.atoms)
        potential_diff = 1.0  # dN
        potential_diff -= self.species_counts[old_species]
        potential_diff -= 0.5 * N * self._concentration_parameters[old_species]
        potential_diff += self.species_counts[new_species]
        potential_diff += 0.5 * N * self._concentration_parameters[new_species]
        potential_diff *= self.variance_parameter
        potential_diff /= N

        potential_diff += self._get_property_change([index], [new_species])

        if self._acceptance_condition(potential_diff):
            self.accepted_trials += 1
            self.update_occupations([index], [new_species])
            self._species_counts[old_species] -= 1
            self._species_counts[new_species] += 1

    def _acceptance_condition(self, potential_diff: float) -> bool:
        """
        Evaluates Metropolis acceptance criterion.

        Parameters
        ----------
        potential_diff
            the change in the thermodynamic potential associated
            with the trial step
        """
        if potential_diff < 0:
            return True
        else:
            return np.exp(-potential_diff / (
                self.boltzmann_constant * self.temperature)) > \
                self._next_random_number()

    @property
    def concentration_parameters(self) -> Dict[int, float]:
        """concentration parameters :math:`\\phi_i`, one for each
        element but their sum must be :math:`-2.0`"""
        return self._concentration_parameters

    @concentration_parameters.setter
    def concentration_parameters(self, concentration_parameters):
        if not isinstance(concentration_parameters, dict):
            raise TypeError('concentration_parameters has the wrong type')
        if abs(sum(concentration_parameters.values()) + 2) > 1e-6:
            raise ValueError('The sum of all concentration_parameters must '
                             'equal -2')

        self._concentration_parameters = {}

        for key, phi in concentration_parameters.items():
            if isinstance(key, str):
                atomic_number = atomic_numbers[key]
                self._concentration_parameters[atomic_number] = phi
            elif isinstance(key, int):
                self._concentration_parameters[key] = phi
        if set(self.configuration._allowed_species) != \
           set(self._concentration_parameters.keys()):
            raise ValueError('concentration_parameters were not set for '
                             'all species')

    @property
    def species_counts(self) -> Dict[int, int]:
        """keeps track of the number of atoms of each species"""
        return self._species_counts

    @species_counts.setter
    def species_counts(self, occupations):
        # Initialize counter to keep track of concentrations efficiently
        species, counts = np.unique(occupations,
                                    return_counts=True)
        self._species_counts = dict(zip(species, counts))

        # There may be species that are not in the input structure
        for species in self.configuration._allowed_species:
            self._species_counts[species] = \
                self._species_counts.get(species, 0)

    def get_ensemble_data(self) -> Dict:
        """
        Returns a dict with the default data of
        the ensemble.

        Here temperature and species counts
        are added to the default data.

        Returns
        -------
        dict : ensemble data key pairs

        """
        data = super().get_ensemble_data()

        # concentration parameters
        for atnum, phi in self.concentration_parameters.items():
            data['phi_{}'.format(chemical_symbols[atnum])] = phi

        # variance parameter
        data['kappa'] = self.variance_parameter

        # free energy derivative
        atnum_1 = min(self.concentration_parameters.keys())
        concentration = self.species_counts[atnum_1] / len(self.atoms)
        data['free_energy_derivative'] = \
            - 2 * self.variance_parameter * concentration - \
            self.variance_parameter * \
            self.concentration_parameters[atnum_1]

        # temperature
        data['temperature'] = self.temperature

        # species counts
        atoms = self.configuration.atoms
        unique, counts = np.unique(atoms.numbers, return_counts=True)
        # TODO: avoid accessing a protected member of a client class
        for atnum in self.configuration._allowed_species:
            data['{}_count'.format(chemical_symbols[atnum])] = 0
        for atnum, count in zip(unique, counts):
            data['{}_count'.format(chemical_symbols[atnum])] = count

        return data

    def update_occupations(self, sites: List[int], species: List[int]):
        """Updates the occupation vector of the configuration being sampled.
        This will change the state of the configuration in both the
        calculator and the configuration manager.

        Parameters
        ----------
        sites
            indices of sites of the configuration to change
        species
            new occupations (species) by atomic number

        Raises
        ------
        ValueError
            if input lists are not of the same length
        """
        super().update_occupations(sites, species)
        self.species_counts = self.configuration.occupations
