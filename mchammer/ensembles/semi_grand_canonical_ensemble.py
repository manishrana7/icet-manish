"""
Definition of the semi-grand canonical ensemble class.
"""

import numpy as np

from ase.data import atomic_numbers, chemical_symbols
from ase.units import kB
from mchammer.ensembles.base_ensemble import BaseEnsemble


class SemiGrandCanonicalEnsemble(BaseEnsemble):
    """Semi-grand canonical (SGC) ensemble.

    The probability density of the SGC ensemble for a
    :math:`m`-component system is

    .. math::
        \\rho_\\text{SGC} = \exp\\Big[ - \\big( E
        + \sum_{i>1}^m \Delta\mu_i c_i \\big) / k_B T \\Big]

    with the *relative* chemical potentials :math:`\\Delta\mu_i =
    \mu_i - \mu_1` and concentrations :math:`c_i = N_i /
    N_\\text{total}`

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

    def __init__(self, atoms=None, calculator=None,
                 name='Semi-grand canonical ensemble',
                 data_container=None, random_seed=None,
                 ensemble_data_write_interval=None, **kwargs):

        super().__init__(atoms=atoms, calculator=calculator, name=name,
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

        if 'chemical_potentials' not in kwargs.keys():
            raise KeyError('Missing required keyword: chemical_potentials')
        else:
            # TODO: check that length of chemical_potentials is correct
            self.chemical_potentials = kwargs['chemical_potentials']

    def do_trial_step(self):
        """ Carry out one Monte Carlo trial step. """
        self.total_trials += 1

        # energy change
        sublattice_index = self.get_random_sublattice_index()
        index, element = \
            self.configuration.get_flip_state(sublattice_index)
        potential_diff = self.get_property_change([index], [element])

        # change in chemical potential
        old_element = self.configuration.occupations[index]
        chemical_potential_diff = \
            self.chemical_potentials[old_element] - \
            self.chemical_potentials[element]
        potential_diff += chemical_potential_diff

        if self._acceptance_condition(potential_diff):
            self.accepted_trials += 1
            self.update_occupations([index], [element])

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
            return np.exp(-potential_diff/(
                self.boltzmann_constant * self.temperature)) > \
                self.next_random_number()

    @property
    def chemical_potentials(self):
        """ dict : chemical potentials :math:`\\mu_i` """
        return self._chemical_potentials

    @chemical_potentials.setter
    def chemical_potentials(self, chemical_potentials):
        # TODO: check that length of chemical_potentials is correct
        self._chemical_potentials = {}
        for key in chemical_potentials.keys():
            if isinstance(key, str):
                element_number = atomic_numbers[key]
                self._chemical_potentials[element_number] =\
                    chemical_potentials[key]
            elif isinstance(key, int):
                self._chemical_potentials[key] =\
                    chemical_potentials[key]

    def get_ensemble_data(self):
        """ 
        Returns a dict with the default data of 
        the ensemble.

        Here the element counts are added to
        the default data.
        """
        default_data = super().get_ensemble_data()

        possible_elements = self.configuration._possible_elements
        atoms = self.configuration.atoms
        unique, counts = np.unique(atoms.numbers, return_counts=True)

        for Z, count in zip(unique, counts):
            str_element = chemical_symbols[Z]
            default_data["{} count".format(str_element)] = count

        # Add the "empty" elements also
        for possible_element in possible_elements:
            if possible_element not in unique:
                str_element = chemical_symbols[possible_element]
                default_data["{} count".format(str_element)] = 0

        return default_data
