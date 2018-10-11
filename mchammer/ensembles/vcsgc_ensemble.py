"""
Definition of the semi-grand canonical ensemble class.
"""

import numpy as np

from ase import Atoms
from ase.data import atomic_numbers, chemical_symbols
from ase.units import kB
from typing import Dict

from .. import DataContainer
from .base_ensemble import BaseEnsemble
from ..calculators.base_calculator import BaseCalculator


class VCSGCEnsemble(BaseEnsemble):
    """Variance-constrained semi-grand canonical (VCSGC) ensemble.
    """

    def __init__(self, atoms: Atoms=None, calculator: BaseCalculator=None,
                 name: str='Semi-grand canonical ensemble',
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

        self._concentration_parameters = concentration_parameters
        self.variance_parameter = variance_parameter

        # Initialize counter to keep track of concentrations efficiently
        species, counts = np.unique(self.configuration.occupations,
                                     return_counts=True)
        self.species_count = dict(zip(species, counts))
        # There may be species that are not in the input structure
        for species in self.configuration._allowed_species:
            self.species_count[species] = self.species_count.get(species, 0)


    def _do_trial_step(self):
        """ Carry out one Monte Carlo trial step. """
        self.total_trials += 1

        # choose swap
        sublattice_index = self.get_random_sublattice_index()
        index, new_species = \
            self.configuration.get_flip_state(sublattice_index)
        old_species = self.configuration.occupations[index]

        # Calculate difference in VCSGC thermodynamic potential
        # Note that this assumes that only one atom was flipped
        potential_diff = 1.0  # dN
        potential_diff -= self.species_count[old_species]
        potential_diff -= 0.5*self._concentration_parameters[old_species]
        potential_diff += self.species_count[new_species]
        potential_diff += 0.5*self._concentration_parameters[new_species]
        potential_diff *= self.variance_parameter
        potential_diff /= len(self.configuration.atoms)

        potential_diff += self._get_property_change([index], [new_species])

        if self._acceptance_condition(potential_diff):
            self.accepted_trials += 1
            self.update_occupations([index], [new_species])
            self.species_count[old_species] -= 1
            self.species_count[new_species] += 1

    def _acceptance_condition(self, potential_diff: float) -> bool:
        """
        Evaluate Metropolis acceptance criterion.

        Parameters
        ----------
        potential_diff : float
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
        """ dict : concentration parameters :math:`\\phi_i` """
        return self._concentration_parameters

    @concentration_parameters.setter
    def concentration_parameters(self, concentration_parameters):
        # TODO: check that length of concentration_parameters is correct
        self._concentration_parameters = {}
        for key, concentration in concentration_parameters.items():
            if isinstance(key, str):
                atomic_number = atomic_numbers[key]
                self._concentration_parameters[atomic_number] =\
                    concentration * len(self.atoms)
            elif isinstance(key, int):
                self._concentration_parameters[key] =\
                    concentration * len(self.atoms)

    def get_ensemble_data(self) -> Dict:
        """
        Returns a dict with the default data of
        the ensemble.

        Here the temperature and species counts
        are added to the default data.

        Returns
        -------
        dict : ensemble data key pairs

        """
        data = super().get_ensemble_data()

        # concentration parameters
        for atnum, phi in self.concentration_parameters.items():
            data['phi_{}'.format(chemical_symbols[atnum])] = phi

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

