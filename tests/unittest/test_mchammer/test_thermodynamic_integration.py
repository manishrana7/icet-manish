from ase.build import bulk
from ase.units import kB
import numpy as np
import unittest

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import ThermodynamicIntegrationEnsemble
from mchammer.free_energy_tools import get_free_energy_thermodynamic_integration
from mchammer.ensembles.thermodynamic_base_ensemble import ThermodynamicBaseEnsemble
from mchammer.ensembles import CanonicalEnsemble

import os
from itertools import permutations
import tempfile

from mchammer.free_energy_tools import _lambda_function_forward


class TestEnsemble(unittest.TestCase):
    """Container for tests of the class functionality."""

    def __init__(self, *args, **kwargs):
        """
        Setup the calculator and calculate the free energy of the fictitious system
        by enumerating all possible occupations
        """
        super(TestEnsemble, self).__init__(*args, **kwargs)
        self.prim = bulk('MoC', crystalstructure='rocksalt',
                         a=4.2)
        self.cs = ClusterSpace(self.prim, cutoffs=[2.2],
                               chemical_symbols=[['Mo', 'V'],
                                                 ['C', 'X', 'O']])
        self.ce = ClusterExpansion(self.cs, [0.0, 0.0, 6.0,  -2.7, 1.2, -0.2])
        self.T = 1200
        self.n_steps = 40000
        self.supercell = self.prim.repeat(2)
        self.natoms = len(self.supercell)
        self.get_all_occupations()
        self.A_enumerated = self.get_free_energy()

    def get_all_occupations(self):
        """ enumerate all possible occupations of the system """
        carbon = np.where(self.supercell.numbers == 6)[0]
        molybdenum = np.where(self.supercell.numbers == 42)[0]
        first_lattice = self.supercell[carbon].numbers
        first_lattice[0] = 8
        first_lattice[1] = 8
        first_lattice[2] = 0
        first_lattice[3] = 0
        second_lattice = self.supercell[molybdenum].numbers
        second_lattice[0] = 23
        second_lattice[1] = 23
        perm1 = list(set(permutations(first_lattice)))
        perm2 = list(set(permutations(second_lattice)))
        self.atoms = []
        for p1 in perm1:
            _atom = self.supercell.copy()
            _atom.numbers[carbon] = p1
            for p2 in perm2:
                atom = _atom.copy()
                atom.numbers[molybdenum] = p2
                self.atoms.append(atom)

    def get_all_energies(self):
        """ get the energy of all possible occupations of the system """
        energies = []
        for atom in self.atoms:
            energies.append(self.ce.predict(atom))
        self.energies = np.array(energies)

    def get_free_energy(self):
        """ calculate the free energy of the system """
        self.get_all_energies()
        beta = 1 / (kB * self.T)
        Z = np.exp(-beta * self.energies).sum()
        pi = np.exp(-beta * self.energies) / Z
        U = (self.energies * pi).sum()
        S = -kB * (pi * np.log(pi)).sum()
        return U - self.T * S

    def get_free_energy_temperature_dependence(self, temperatures):
        free_energies = []
        for T in temperatures:
            beta = 1 / (kB * T)
            Z = np.exp(-beta * self.energies).sum()
            pi = np.exp(-beta * self.energies) / Z
            U = (self.energies * pi).sum()
            S = -kB * (pi * np.log(pi)).sum()
            free_energies.append(U - T * S)
        return free_energies

    def setUp(self):
        """ setup the calculator """
        self.calculator = ClusterExpansionCalculator(self.atoms[0], self.ce,
                                                     scaling=1.0)
        ensemble = CanonicalEnsemble(
                self.atoms[0], calculator=self.calculator,
                temperature=self.T,
                dc_filename=None)
        ensemble.run(60000)
        self.temperature_structure = ensemble.structure

    def test_run_forward(self):
        """
        run a thermodynamic integration simulation and calculate the free energy
        this is the compared with the direct method
        """
        ensemble = ThermodynamicIntegrationEnsemble(
            structure=self.atoms[0], calculator=self.calculator,
            temperature=self.T,
            n_steps=self.n_steps,
            forward=True,
            dc_filename=None)
        ensemble.run()
        dc = ensemble.data_container
        (_, A) = get_free_energy_thermodynamic_integration(dc, self.cs, True, self.T)
        assert np.allclose(A[0], self.A_enumerated, atol=5e-3)

    def test_run_backward(self):
        """
        run a thermodynamic integration simulation and calculate the free energy
        this is the compared with the direct method
        """
        ensemble = ThermodynamicIntegrationEnsemble(
            structure=self.temperature_structure, calculator=self.calculator,
            temperature=self.T,
            n_steps=self.n_steps,
            forward=False,
            dc_filename=None)
        ensemble.run()
        dc = ensemble.data_container
        (temperatures, A) = get_free_energy_thermodynamic_integration(dc, self.cs, False, self.T)
        A_enumerated = self.get_free_energy_temperature_dependence(temperatures)
        assert np.allclose(A, A_enumerated, atol=5e-3)

    def test_run_with_temperature_dependence_foward(self):
        """
        run a thermodynamic integration simulation and calculate the free energy
        this is the compared with the direct method
        """
        ensemble = ThermodynamicIntegrationEnsemble(
            structure=self.atoms[0], calculator=self.calculator,
            temperature=self.T,
            n_steps=self.n_steps,
            forward=True,
            dc_filename=None)
        ensemble.run()
        dc = ensemble.data_container
        (temperatures, A) = get_free_energy_thermodynamic_integration(dc, self.cs,
                                                                      True, 1000)
        A_enumerated = self.get_free_energy_temperature_dependence(temperatures)
        assert np.allclose(A, A_enumerated, atol=5e-3)

    def test_run_with_temperature_dependence_backward(self):
        """
        run a thermodynamic integration simulation and calculate the free energy
        this is the compared with the direct method
        """
        ensemble = ThermodynamicIntegrationEnsemble(
            structure=self.temperature_structure, calculator=self.calculator,
            temperature=self.T,
            n_steps=self.n_steps,
            forward=False,
            dc_filename=None)
        ensemble.run()
        dc = ensemble.data_container
        (temperatures, A) = get_free_energy_thermodynamic_integration(dc, self.cs,
                                                                      False, 1000)
        A_enumerated = self.get_free_energy_temperature_dependence(temperatures)
        assert np.allclose(A, A_enumerated, atol=5e-3)

    def test_restart_run(self):
        """
        run a thermodynamic integration simulation and calculate the free energy
        this is the compared with the direct method.
        However, this time the simulation is restarted.
        """
        tmpfolder = tempfile._get_default_tempdir()
        tmpname = next(tempfile._get_candidate_names())
        fname = f'{tmpfolder}/{tmpname}'
        ensemble = ThermodynamicIntegrationEnsemble(
            structure=self.atoms[0], calculator=self.calculator,
            temperature=self.T,
            n_steps=self.n_steps,
            forward=True,
            dc_filename=fname)

        ensemble._lambda_function = _lambda_function_forward
        ThermodynamicBaseEnsemble.run(ensemble, 100)

        new_ensemble = ThermodynamicIntegrationEnsemble(
            structure=self.atoms[0], calculator=self.calculator,
            temperature=self.T,
            n_steps=self.n_steps,
            forward=True,
            dc_filename=fname)
        new_ensemble.run()
        dc = new_ensemble.data_container
        (_, A) = get_free_energy_thermodynamic_integration(dc, self.cs, True, self.T)
        os.remove(fname)
        assert np.allclose(A[0], self.A_enumerated, atol=5e-3)
