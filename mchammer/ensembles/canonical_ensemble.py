"""Definition of the canonical ensemble class."""

import numpy as np
from ase import Atoms
from ase.units import kB
from typing import Dict

from .. import DataContainer
from .base_ensemble import BaseEnsemble
from ..calculators.base_calculator import BaseCalculator


class CanonicalEnsemble(BaseEnsemble):
    """Canonical Ensemble.

    Instances of this class allow one to simulate systems in the
    canonical ensemble (:math:`N_iVT`), i.e. at constant temperature
    (:math:`T`), number of species (:math:`N_i`), and volume
    (:math:`V`).

    Attributes
    -----------
    temperature : float
        temperature :math:`T` in appropriate units [commonly Kelvin]
    boltzmann_constant : float
        Boltzmann constant :math:`k_B` in appropriate
        units, i.e. units that are consistent
        with the underlying cluster expansion
        and the temperature units [default: eV/K]

    """

    def __init__(self, atoms: Atoms=None, calculator: BaseCalculator=None,
                 name: str='Canonical ensemble',
                 data_container: DataContainer=None, random_seed: int=None,
                 ensemble_data_write_interval: int=None, **kwargs):
        super().__init__(
            atoms=atoms, calculator=calculator, name=name,
            data_container=data_container,
            random_seed=random_seed,
            ensemble_data_write_interval=ensemble_data_write_interval)
        if 'temperature' not in kwargs.keys():
            raise KeyError('Missing required keyword: temperature')
        else:
            self.temperature = kwargs['temperature']

        if 'boltzmann_constant' in kwargs.keys():
            self.boltzmann_constant = kwargs['boltzmann_constant']
        else:
            self.boltzmann_constant = kB

    def do_trial_step(self):
        """ Carries out one Monte Carlo trial step. """
        self.total_trials += 1

        sublattice_index = self.get_random_sublattice_index()
        sites, species = \
            self.configuration.get_swapped_state(sublattice_index)

        potential_diff = self.get_property_change(sites, species)

        if self._acceptance_condition(potential_diff):
            self.accepted_trials += 1
            self.update_occupations(sites, species)

    def _acceptance_condition(self, potential_diff: float) -> bool:
        """
        Evaluates Metropolis acceptance criterion.

        Parameters
        ----------
        potential_diff
            change in the thermodynamic potential associated
            with the trial step
        """
        if potential_diff < 0:
            return True
        else:
            return np.exp(-potential_diff/(
                self.boltzmann_constant * self.temperature)) > \
                self.next_random_number()

    def get_ensemble_data(self) -> Dict:
        """Returns the data associated with the ensemble. For the SGC
        ensemble this specifically includes the temperature and the
        species counts.
        """
        data = super().get_ensemble_data()
        data['temperature'] = self.temperature
        return data
