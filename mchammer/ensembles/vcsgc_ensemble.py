"""
Definition of the semi-grand canonical ensemble class.
"""

import numpy as np

from ase.data import atomic_numbers, chemical_symbols
from ase.units import kB
from mchammer.ensembles.base_ensemble import BaseEnsemble

from typing import Dict


class VCSGCEnsemble(BaseEnsemble):
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
                 name='VCSGC ensemble',
                 data_container=None, random_seed=None,
                 ensemble_data_write_interval=None, **kwargs):

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

        if 'concentration_parameters' not in kwargs.keys():
            raise KeyError(
                'Missing required keyword: concentration_parameters')
        else:
            # TODO: check that length of concentration_parameter is correct
            self.concentration_parameters = kwargs['concentration_parameters']

        if 'variance_parameter' not in kwargs.keys():
            raise KeyError('Missing required keyword: variance_parameter')
        else:
            self.variance_parameter = kwargs['variance_parameter']

        # Initialize counter to keep track of concentrations efficiently
        elements, counts = np.unique(self.configuration.occupations,
                                     return_counts=True)
        self.element_count = dict(zip(elements, counts))
        # There may be elements that are not in the input structure
        for element in self.configuration._possible_elements:
            self.element_count[element] = self.element_count.get(element, 0)


    def do_trial_step(self):
        """ Carry out one Monte Carlo trial step. """
        self.total_trials += 1

        # choose swap
        sublattice_index = self.get_random_sublattice_index()
        index, new_element = \
            self.configuration.get_flip_state(sublattice_index)
        old_element = self.configuration.occupations[index]

        # Calculate difference in VCSGC thermodynamic potential
        # Note that this assumes that only one atom was flipped
        potential_diff = 1.0  # dN
        potential_diff += self.element_count[old_element]
        potential_diff -= self.concentration_parameters[old_element]
        potential_diff -= self.element_count[new_element]
        potential_diff += self.concentration_parameters[new_element]
        potential_diff *= 2 * self.variance_parameter
        potential_diff /= len(self.configuration.atoms)

        potential_diff += self.get_property_change([index], [new_element])

        if self._acceptance_condition(potential_diff):
            self.accepted_trials += 1
            self.update_occupations([index], [new_element])
            self.element_count[old_element] -= 1
            self.element_count[new_element] += 1

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
                self.next_random_number()

    @property
    def concentration_parameters(self):
        """ dict : concentration parameters :math:`\\phi_i` """
        return self._concentration_parameters

    @concentration_parameters.setter
    def concentration_parameters(self, concentration_parameters):
        # TODO: check that length of concentration_parameters is correct
        self._concentration_parameters = {}
        for key, concentration in concentration_parameters.items():
            if isinstance(key, str):
                element_number = atomic_numbers[key]
                self._concentration_parameters[element_number] =\
                    len(self.configuration.atoms) * concentration
            elif isinstance(key, int):
                self._concentration_parameters[key] =\
                    len(self.configuration.atoms) * concentration

    def get_ensemble_data(self) -> Dict:
        """
        Returns a dict with the default data of
        the ensemble.

        Here the temperature and element counts
        are added to the default data.

        Returns
        -------
        dict : ensemble data key pairs

        """
        data = super().get_ensemble_data()
        data['temperature'] = self.temperature

        atoms = self.configuration.atoms
        unique, counts = np.unique(atoms.numbers, return_counts=True)

        for Z, count in zip(unique, counts):
            str_element = chemical_symbols[Z]
            data["{} count".format(str_element)] = count

        # Add the "empty" elements also
        for possible_element in self.configuration._possible_elements:
            if possible_element not in unique:
                str_element = chemical_symbols[possible_element]
                data["{} count".format(str_element)] = 0

        return data

    def test(self):
        print(self.configuration)
