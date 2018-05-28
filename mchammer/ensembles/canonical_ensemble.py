"""Definition of the canonical ensemble class."""

from mchammer.ensembles.base_ensemble import BaseEnsemble
import numpy as np
from ase.units import kB


class CanonicalEnsemble(BaseEnsemble):
    """
    Canonical Ensemble.

    Attributes
    -----------
    temperature : float
        temperature in Kelvin.
    boltzmann_constant : float
        Boltzmann constant in appropriate
        units, i.e. units that are consistent
        with the underlying cluster expansion
        and the temperature [default: eV/K]
    """

    def __init__(self, atoms=None, calculator=None, name='Canonical Ensemble',
                 data_container=None, random_seed=None,ensemble_data_write_interval=None, **kwargs):

        super().__init__(atoms=atoms, calculator=calculator, name=name,
                         data_container=data_container,
                         random_seed=random_seed,
                         ensemble_data_write_interval=ensemble_data_write_interval)
        if 'temperature' not in kwargs.keys():
            raise KeyError("Temperature needs to be set in canonical ensemble")
        else:
            self.temperature = kwargs['temperature']
        if 'boltzmann_constant' in kwargs.keys():
            self.boltzmann_constant = kwargs['boltzmann_constant']
        else:
            self.boltzmann_constant = kB

    def do_trial_step(self):
        """Do a trial step."""
        self.total_trials += 1

        sublattice_index = self.get_random_sublattice_index()
        indices, elements = \
            self.configuration.get_swapped_state(sublattice_index)

        energy_diff = self.get_property_change(indices, elements)

        if self._acceptance_condition(energy_diff):
            self.accepted_trials += 1
            self.update_occupations(indices, elements)

    def _acceptance_condition(self, energy_diff: float) -> bool:
        """
        Evaluate Metropolis acceptance criterion.

        Parameters
        ----------
        energy_diff : float
            the energy difference associated with
            this trial step
        """
        if energy_diff < 0:
            return True
        else:
            return np.exp(-energy_diff/(
                self.boltzmann_constant * self.temperature)) > \
                self.next_random_number()
