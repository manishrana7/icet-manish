from typing import List, Union, Dict

from icet import ClusterSpace
from ase import Atoms
from ase.data import chemical_symbols as chemical_symbols_ase
import numpy as np
import copy

from icet.core.structure import Structure
from _icet import ClusterSpace as _ClusterSpace


class VariableBasisClusterSpace(ClusterSpace):

    def __init__(self,
                 structure: Atoms,
                 cutoffs: List[float],
                 chemical_symbols: Union[List[str], List[List[str]]],
                 nparameters_per_orbit: List[float],
                 symprec: float = 1e-5,
                 position_tolerance: float = None):
        super().__init__(structure=structure,
                         cutoffs=cutoffs,
                         chemical_symbols=chemical_symbols,
                         point_function_form='hyperbolic',
                         symprec=symprec,
                         position_tolerance=position_tolerance)

        # Ensure that there is only one active sublattice
        active_sublattices = self.get_sublattices(structure).active_sublattices
        if len(active_sublattices) != 1:
            raise NotImplementedError('VariableBasisClusterSpace can currently only be used with'
                                      ' one active sublattice')
        self._active_chemical_symbols = tuple(active_sublattices[0].chemical_symbols)

        # Remove the singlet (it is always identically zero with the variable basis)
        self._remove_orbits([0])

        # Ensure that the number of parameters per orbit matches
        # the length of the (standard) cluster space
        if len(nparameters_per_orbit) != len(self._orbit_list) + 1:  # +1 because zerolet
            raise ValueError(f"The length of nparameters_per_orbit ({len(nparameters_per_orbit)})"
                             " must be equal to length of the (standard) cluster space,"
                             f" excluding the singlet ({len(self._orbit_list) + 1}).")
        self._nparameters_per_orbit = nparameters_per_orbit

        # Extract a list of number of bodies per orbit
        nbodies_per_orbit = []
        for orbit in self.orbit_data:
            nbodies_per_orbit.append(orbit['order'])
        assert len(nbodies_per_orbit) == len(self._nparameters_per_orbit)
        self._nbodies_per_orbit = nbodies_per_orbit

        # Decide how to determine which is atom type A when calculating
        # concentration as N_A / N. This decision needs to be consistent
        # with how icet decides spin variables internally.
        atomic_number = min(active_sublattices[0].atomic_numbers)
        self._concentration_symbol = chemical_symbols_ase[atomic_number]

    def get_cluster_vector(self, structure: Atoms) -> np.ndarray:
        """
        Returns the cluster vector for a structure.

        TODO: Consider speeding up this in some smart way.

        Parameters
        ----------
        structure
            atomic configuration

        Returns
        -------
        the cluster vector
        """
        if not isinstance(structure, Atoms):
            raise TypeError('Input structure must be an ASE Atoms object')

        x = 1 - 2 * self.get_concentration(structure)
        try:
            cv = _ClusterSpace.get_cluster_vector(
                self,
                structure=Structure.from_atoms(structure),
                fractional_position_tolerance=self.fractional_position_tolerance,
                point_function_parameter=x
                )
        except Exception as e:
            self.assert_structure_compatibility(structure)
            raise(e)
        
        # Now extend the cluster vector with the desired
        # number of parameters per element
        assert len(cv) == len(self._nparameters_per_orbit)
        extended_cv = []
        for i, corr in enumerate(cv):
            for m in range(self._nparameters_per_orbit[i]):
                prefactor = get_cluster_correlation_prefactor(x, m, self._nbodies_per_orbit[i])
                extended_cv.append(prefactor * corr)
        return extended_cv

    def get_concentration(self, structure: Atoms):
        """
        Get concentration for a structure.

        TODO: The present implementation is clunky and probably very slow.
        """
        active_sublattices = self.get_sublattices(structure).active_sublattices
        assert len(active_sublattices) == 1
        active_sublattice = active_sublattices[0]
        symbols = np.array(structure.get_chemical_symbols())[active_sublattice.indices]
        return list(symbols).count(self.concentration_symbol) / len(symbols)

    @property
    def concentration_symbol(self):
        return self._concentration_symbol

    @property
    def active_chemical_symbols(self):
        return copy.copy(self._active_chemical_symbols)

    @property
    def nparameters_per_orbit(self):
        return copy.copy(self._nparameters_per_orbit)

    @property
    def total_number_of_parameters(self):
        return sum(self._nparameters_per_orbit)

    def _get_cluster_space_parameters_dict(self,
                                           include_input_structure: bool = False,
                                           include_pruning_history: bool = False) -> Dict:
        """
        Returns a dict with parameters that specify this cluster space
        (used for writing cluster space to file).
        """
        items = super()._get_cluster_space_parameters_dict(
                include_input_structure=include_input_structure,
                include_pruning_history=include_pruning_history)
        items['nparameters_per_orbit'] = self._nparameters_per_orbit
        return items


def get_cluster_correlation_prefactor(x: float, m: int, n: int):
    return (1 - x**2)**(n / 2) * x**m
