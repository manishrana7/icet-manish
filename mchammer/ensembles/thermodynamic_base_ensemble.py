import numpy as np

from ase import Atoms
from ase.units import kB
from typing import Dict

from .. import DataContainer
from .base_ensemble import BaseEnsemble
from ..calculators.base_calculator import BaseCalculator
from abc import abstractproperty

class ThermodynamicBaseEnsemble(BaseEnsemble):
    """
    Parameters
    ----------
    atoms : :class:`Atoms <ase.Atoms>`
        atomic configuration to be used in the Monte Carlo simulation;
        also defines the initial occupation vector
    calculator : :class:`BaseCalculator <mchammer.calculators.ClusterExpansionCalculator>`
        calculator to be used for calculating the potential changes
        that enter the evaluation of the Metropolis criterion
    T_start : float
        temperature from which the annealing is started
    T_stop : float
        final temperature for annealing
    n_steps : int
        number of steps to take in the annealing simulation
    cooling_function : str/function
        to use the predefined cooling functions provide a string
        `linear` or `exponential`, otherwise provide a function
    boltzmann_constant : float
        Boltzmann constant :math:`k_B` in appropriate
        units, i.e. units that are consistent
        with the underlying cluster expansion
        and the temperature units [default: eV/K]
    user_tag : str
        human-readable tag for ensemble [default: None]
    data_container : str
        name of file the data container associated with the ensemble
        will be written to; if the file exists it will be read, the
        data container will be appended, and the file will be
        updated/overwritten
    random_seed : int
        seed for the random number generator used in the Monte Carlo
        simulation
    ensemble_data_write_interval : int
        interval at which data is written to the data container; this
        includes for example the current value of the calculator
        (i.e. usually the energy) as well as ensembles specific fields
        such as temperature or the number of atoms of different species
    data_container_write_period : float
        period in units of seconds at which the data container is
        written to file; writing periodically to file provides both
        a way to examine the progress of the simulation and to back up
        the data [default: np.inf]
    trajectory_write_interval : int
        interval at which the current occupation vector of the atomic
        configuration is written to the data container.
    """

    def __init__(self, atoms: Atoms, calculator: BaseCalculator,
                 user_tag: str = None,
                boltzmann_constant: float = kB,
                data_container: DataContainer = None, random_seed: int = None,
                data_container_write_period: float = np.inf,
                ensemble_data_write_interval: int = None,
                trajectory_write_interval: int = None) -> None:
        
        
        self._boltzmann_constant = boltzmann_constant

        super().__init__(
            atoms=atoms, calculator=calculator, user_tag=user_tag,
            data_container=data_container,
            random_seed=random_seed,
            data_container_write_period=data_container_write_period,
            ensemble_data_write_interval=ensemble_data_write_interval,
            trajectory_write_interval=trajectory_write_interval)



    @abstractproperty
    @property
    def temperature(self) -> float:
        pass

    @property
    def boltzmann_constant(self) -> float:
        """ Boltzmann constant :math:`k_B` (see parameters section above) """
        return self._boltzmann_constant


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
        elif abs(self.temperature) < 1e-6:  # temperature is numerically zero
            return False
        else:
            p = np.exp(-potential_diff / (self.boltzmann_constant * self.temperature))
            return p > self._next_random_number()

    def do_canonical_swap(self, sublattice_index=None):
        """ Carries out one Monte Carlo trial step. """


        sites, species = self.configuration.get_swapped_state(sublattice_index)

        potential_diff = self._get_property_change(sites, species)

        if self._acceptance_condition(potential_diff):
            self._accepted_trials += 1
            self.update_occupations(sites, species)


    def do_sgc_flip(self, chemical_potentials, sublattice_index=None):
        """ Carries out one Monte Carlo trial step. """

        if sublattice_index == None:
            sublattice_index = self.get_random_sublattice_index()

        index, species = \
            self.configuration.get_flip_state(sublattice_index)
        potential_diff = self._get_property_change([index], [species])

        # change in chemical potential
        old_species = self.configuration.occupations[index]
        chemical_potential_diff = \
            chemical_potentials[old_species] - \
            chemical_potentials[species]
        potential_diff += chemical_potential_diff

        if self._acceptance_condition(potential_diff):
            self._accepted_trials += 1
            self.update_occupations([index], [species])
