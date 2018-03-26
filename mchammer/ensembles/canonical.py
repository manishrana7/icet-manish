"""Definition of the canonical ensemble class."""

from mchammer.ensembles.base_ensemble import BaseEnsemble
import numpy as np
from ase.units import kB

class Canonical(BaseEnsemble):
    """
    Canonical Ensemble.
    
    Attributes:
    temperature : float
        temperature in kelvin.
    """

    def __init__(self, atoms, temperature, calculator, name='Canonical Ensemble',
                 data_container=None, random_seed=None):

        super().__init__(atoms=atoms, calculator=calculator, name=name,
                         data_container=data_container,
                         random_seed=random_seed)

        self.temperature = temperature


    def do_trial_move(self):
        """
        Do a trial move.
        """
        
        index1, index2 = self._get_swap_indices()

        energy_before = self.calculate_property([index1, index2])
        self._canonical_swap(index1, index2)

        energy_after = self.calculate_property([index1, index2])
        energy_diff = energy_after - energy_before

        self.total_trials +=1
        if energy_diff < 0 or self._acceptance_condition(energy_diff):
            self.accepted_trials +=1
        
            
    def _acceptance_condition(self, energy_diff):
        """
        Check if this energy diff is accepted.
        """

        return np.exp(-energy_diff/(kB * self.temperature)) > self.next_random_number()

    def _run(self, number_of_trial_moves):
        """
        Private run method
        """
        pass

    def _get_swap_indices(self):
        """
        Get indices for swapping.
        """

        return 0,1

    def _canonical_swap(self, index1, index2):
        """
        Swap atoms on indices  index1 and index2.

        Parameters:
        index1 : int
            first integer to swap.
        index2 : int
            second integer to swap.
        """

        symbol1 = self.structure[index1].symbol
        symbol2 = self.structure[index2].symbol

        self.structure[index1].symbol = symbol2
        self.structure[index2].symbol = symbol1