from ase.units import kB
import numpy as np
import unittest
from ase import Atoms
from ase.build import bulk

from icet import ClusterSpace
from mchammer.free_energy_tools import (get_free_energy_thermodynamic_integration,
                                        get_free_energy_temperature_integration)
from mchammer import DataContainer
from mchammer.free_energy_tools import (_lognpermutations, _npermutations)

from scipy.special import binom
from math import factorial


class TestEnsemble(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        super(TestEnsemble, self).__init__(*args, **kwargs)
        self.prim_single = Atoms('Au', positions=[[0, 0, 0]], cell=[1, 1, 10], pbc=True)
        self.cs_single = ClusterSpace(self.prim_single, cutoffs=[1.1],
                                      chemical_symbols=['Ag', 'Au'])

        self.prim_double = bulk('MoC', crystalstructure='rocksalt',
                                a=4.2)
        self.cs_double = ClusterSpace(self.prim_double, cutoffs=[2.2],
                                      chemical_symbols=[['Mo', 'V'],
                                                        ['C', 'X', 'O']])

    def test_thermodynamic_integration_single_sublattice(self):
        """
        Fictive simple test for the calculation of the free energy from
        a themodynamic integration simulation when the potential is constantly 0.3
        such that we know the analytical value of the integration.
        There is two Ag to be placed on n_atoms such that the boltzmann entropy
        is n_atoms choose 2.
        """

        structure = self.prim_single.repeat(5)
        n_Ag = 2
        for n in range(n_Ag):
            structure[n].symbol = 'Ag'
        temperature = 100
        potential = 0.3
        dc = DataContainer(
                structure=structure,
                ensemble_parameters={'n_atoms': len(structure),
                                     'temperature': temperature})
        for i in range(1000):
            dc.append(i, {'potential': potential,
                          'occupations': structure.numbers})
        (_, A) = get_free_energy_thermodynamic_integration(dc, self.cs_single,
                                                           True,
                                                           temperature)
        s = -kB * np.log(binom(len(structure), n_Ag))
        A_direct = temperature * s + potential
        assert np.allclose(A_direct, A[0], atol=1e-5)

    def test_thermodynamic_integration_multiple_sublattice(self):
        """
        Fictive test for the calculation of the free energy from
        a themodynamic integration simulation when the potential is constantly 0.3
        such that we know the analytical value of the integration.

        The first sublattice contains 3 different atomic types C,O,X,
        whereas the second sublattice contains 2 different atomic types, Mo,V

        This test checks that the number of permutations for the two sublattices
        is correctly calculated
        """

        structure = self.prim_double.repeat(5)
        carbon = np.where(structure.numbers == 6)[0]
        molybdenum = np.where(structure.numbers == 42)[0]
        n_vacancies = 5
        for n in range(n_vacancies):
            structure[carbon[n]].symbol = 'X'
        n_oxygen = 3
        for n in range(n_vacancies, n_vacancies + n_oxygen):
            structure[carbon[n]].symbol = 'O'
        n_vanadium = 6
        for n in range(n_vanadium):
            structure[molybdenum[n]].symbol = 'V'

        temperature = 100
        potential = 0.3
        dc = DataContainer(
                structure=structure,
                ensemble_parameters={'n_atoms': len(structure),
                                     'temperature': temperature})
        for i in range(1000):
            dc.append(i, {'potential': potential,
                          'occupations': structure.numbers})
        (_, A) = get_free_energy_thermodynamic_integration(dc, self.cs_double,
                                                           True,
                                                           temperature)
        s = kB * (np.log(factorial(len(molybdenum)) /
                         (factorial(len(molybdenum) - n_vanadium) * factorial(n_vanadium))) +
                  np.log(factorial(len(carbon)) /
                         (factorial(len(carbon) - n_vacancies - n_oxygen)
                          * factorial(n_vacancies) * factorial(n_oxygen))))
        A_direct = potential - temperature * s
        assert np.allclose(A_direct, A[0], atol=1e-5)

    def test_thermodynamic_integration_multiple_sublattice_only_one_active(self):
        """
        Fictive test for the calculation of the free energy from
        a themodynamic integration simulation when the potential is constantly 0.3
        such that we know the analytical value of the integration.

        The first sublattice contains 3 different atomic types C,O,X,
        whereas the second sublattice contains 2 different atomic types, Mo,V
        However, the second lattice has sublattice_probability of 0.

        This test checks that the number of permutations for the two sublattices
        is correctly calculated
        """

        structure = self.prim_double.repeat(5)
        carbon = np.where(structure.numbers == 6)[0]
        molybdenum = np.where(structure.numbers == 42)[0]
        n_vacancies = 5
        for n in range(n_vacancies):
            structure[carbon[n]].symbol = 'X'
        n_oxygen = 3
        for n in range(n_vacancies, n_vacancies + n_oxygen):
            structure[carbon[n]].symbol = 'O'
        n_vanadium = 6
        for n in range(n_vanadium):
            structure[molybdenum[n]].symbol = 'V'

        temperature = 100
        potential = 0.3
        dc = DataContainer(
                structure=structure,
                ensemble_parameters={'n_atoms': len(structure),
                                     'temperature': temperature})
        for i in range(1000):
            dc.append(i, {'potential': potential,
                          'occupations': structure.numbers})
        (_, A) = get_free_energy_thermodynamic_integration(dc, self.cs_double,
                                                           True,
                                                           temperature,
                                                           sublattice_probabilities=[1, 0])
        s = kB * (np.log(factorial(len(carbon)) /
                         (factorial(len(carbon) - n_vacancies - n_oxygen)
                          * factorial(n_vacancies) * factorial(n_oxygen))))
        A_direct = potential - temperature * s
        assert np.allclose(A_direct, A[0], atol=1e-5)

    def test_thermodynamic_integration_multiple_sublattice_only_one_active_temperature(self):
        """
        Fictive test for the calculation of the free energy from
        a themodynamic integration simulation when the potential is constantly 0.3
        such that we know the analytical value of the integration.

        The first sublattice contains 3 different atomic types C,O,X,
        whereas the second sublattice contains 2 different atomic types, Mo,V
        However, the second lattice has sublattice_probability of 0.

        This test checks that the number of permutations for the two sublattices
        is correctly calculated
        """

        structure = self.prim_double.repeat(5)
        carbon = np.where(structure.numbers == 6)[0]
        molybdenum = np.where(structure.numbers == 42)[0]
        n_vacancies = 5
        for n in range(n_vacancies):
            structure[carbon[n]].symbol = 'X'
        n_oxygen = 3
        for n in range(n_vacancies, n_vacancies + n_oxygen):
            structure[carbon[n]].symbol = 'O'
        n_vanadium = 6
        for n in range(n_vanadium):
            structure[molybdenum[n]].symbol = 'V'

        temperature = 100
        potential = 0.3
        dc = DataContainer(
                structure=structure,
                ensemble_parameters={'n_atoms': len(structure),
                                     'temperature': temperature})
        for i in range(1000):
            dc.append(i, {'potential': potential,
                          'occupations': structure.numbers})
        (_, A) = get_free_energy_thermodynamic_integration(dc, self.cs_double,
                                                           True,
                                                           temperature,
                                                           sublattice_probabilities=[1, 0])
        s = kB * (np.log(factorial(len(carbon)) /
                         (factorial(len(carbon) - n_vacancies - n_oxygen)
                          * factorial(n_vacancies) * factorial(n_oxygen))))
        A_direct = potential - temperature * s
        assert np.allclose(A_direct, A[0], atol=1e-5)

    def test_n_permutations_consistency(self):
        """
        This test just checks the consistency of
        _lognpermutations and _npermutations
        _lognpermutations uses stirlings approximation if we get an
        overflow in _npermtations

        The approximation is not 100% accurate so we need to allow for some
        error here
        """
        structure = self.prim_double.repeat(13)
        carbon = np.where(structure.numbers == 6)[0]
        molybdenum = np.where(structure.numbers == 42)[0]
        n_vacancies = 203
        for n in range(n_vacancies):
            structure[carbon[n]].symbol = 'X'
        n_vanadium = 200
        for n in range(n_vanadium):
            structure[molybdenum[n]].symbol = 'V'

        vacancies = np.where(structure.numbers == 0)[0]
        carbon = np.where(structure.numbers == 6)[0]
        molybdenum = np.where(structure.numbers == 42)[0]
        vanadium = np.where(structure.numbers == 23)[0]

        first_lattice = np.append(structure.numbers[molybdenum],
                                  structure.numbers[vanadium])
        second_lattice = np.append(structure.numbers[vacancies],
                                   structure.numbers[carbon])
        logperms_no_approx = (np.log(_npermutations(first_lattice)) +
                              np.log(_npermutations(second_lattice)))
        logperms_approx = (_lognpermutations(first_lattice) +
                           _lognpermutations(second_lattice))
        assert np.allclose(logperms_no_approx, logperms_approx, atol=3)

    def test_temperature_integration_forward(self):
        """
        Fictive simple test for the calculation of the free energy from
        a temperature integration when the potential given as T^2
        such that we know the analytical value of the integration.

        We skip 1000 temperatures since the numerical stability of the
        integration when the number of sample points is small is not great.
        """
        temperatures = np.linspace(50, 10, 10000)
        A_reference = 0
        structure = self.prim_single.repeat(5)
        dc = DataContainer(structure=structure,
                           ensemble_parameters={'n_atoms': len(structure)})
        dc._observables.add('potential')
        dc._observables.add('mctrial')
        dc._observables.add('temperature')
        for i in range(temperatures.size):
            dc._data_list.append({'mctrial': i, 'potential': temperatures[i]**2,
                                  'temperature': temperatures[i]})

        (_, A) = get_free_energy_temperature_integration(dc,
                                                         self.cs_single,
                                                         forward=True,
                                                         temperature_reference=temperatures[0],
                                                         free_energy_reference=A_reference,
                                                         max_temperature=temperatures[-1000],
                                                         )
        integration = (temperatures - temperatures[0])
        A_direct = (A_reference / temperatures[0] - integration) * temperatures

        assert np.allclose(A_direct[-1000:], A, atol=1e-5)

    def test_temperature_integration_backward(self):
        """
        Fictive simple test for the calculation of the free energy from
        a temperature integration when the potential given as T^2
        such that we know the analytical value of the integration.

        We skip 1000 temperatures since the numerical stability of the
        integration when the number of sample points is small is not great.
        """
        temperatures = np.linspace(10, 50, 10000)
        A_reference = 0
        structure = self.prim_single.repeat(5)
        dc = DataContainer(structure=structure,
                           ensemble_parameters={'n_atoms': len(structure)})
        dc._observables.add('potential')
        dc._observables.add('mctrial')
        dc._observables.add('temperature')
        for i in range(temperatures.size):
            dc._data_list.append({'mctrial': i, 'potential': temperatures[i]**2,
                                  'temperature': temperatures[i]})

        (_, A) = get_free_energy_temperature_integration(dc,
                                                         self.cs_single,
                                                         forward=False,
                                                         temperature_reference=temperatures[-1],
                                                         free_energy_reference=A_reference,
                                                         max_temperature=temperatures[1000],
                                                         )
        integration = (temperatures - temperatures[-1])
        A_direct = (A_reference / temperatures[-1] - integration) * temperatures

        assert np.allclose(A_direct[0:1001][::-1], A, atol=1e-5)
