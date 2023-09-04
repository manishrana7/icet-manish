from ase import Atoms
from ase.units import kB
import numpy as np
import matplotlib.pyplot as plt

from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators import ClusterExpansionCalculator
from mchammer.ensembles import ThermodynamicIntegrationEnsemble
from mchammer.ensembles import CanonicalAnnealing
from mchammer.ensembles import CanonicalEnsemble
from mchammer.free_energy_tools import (get_free_energy_thermodynamic_integration,
                                        get_free_energy_temperature_integration)
from sympy.utilities.iterables import multiset_permutations


def get_all_structures(struct):
    """
    Generates all possible lattice occupations
    """
    symbols = struct.get_chemical_symbols()
    for i in range(8):
        symbols[i] = 'Ag'
    permutation = multiset_permutations(symbols)
    atoms = []
    for perm in permutation:
        atom = struct.copy()
        atom.symbols = perm
        atoms.append(atom)
    return atoms


def get_free_energy_enumeration(temperatures, atoms, calc, boltzmann_constant=kB):
    """
    Calculate the free energy by brute force
    """
    energies = []
    for atom in atoms:
        energies.append(calc.calculate_total(occupations=atom.numbers))
    energies = np.array(energies)

    free_energy = []
    for T in temperatures:
        beta = 1 / (boltzmann_constant * T)
        Z = np.exp(-beta * energies).sum()
        free_energy.append(-boltzmann_constant * T * np.log(Z))
    return np.array(free_energy)


##########################################################################################
# Setup
##########################################################################################
# Prepare cluster expansion
prim = Atoms('Au', positions=[[0, 0, 0]], cell=[1, 1, 10], pbc=True)
cs = ClusterSpace(prim, cutoffs=[1.01], chemical_symbols=['Ag', 'Au'])
ce = ClusterExpansion(cs, [0, 0, 2])

# Prepare initial configuration and calculator
supercell = prim.repeat((4, 4, 1))
calc = ClusterExpansionCalculator(supercell, ce)

start_configuration = supercell.copy()
for i in range(8):
    start_configuration[i].symbol = 'Ag'

n_integration_steps = 400000
n_equilibration_steps = 1000
temperature_max = 30
temperature_min = 0.2
temperature_max_plot_limit = 4

##########################################################################################
# Thermodynamic integration
##########################################################################################
# Prepare for the TI by running the canonical ensemble at a high temperature.
mc = CanonicalEnsemble(
        structure=start_configuration,
        calculator=calc,
        temperature=temperature_max,
        boltzmann_constant=1,
        trajectory_write_interval=None,
        ensemble_data_write_interval=200)
mc.run(n_equilibration_steps)

# Run the TI from the disordered to the ordered system.
mc = ThermodynamicIntegrationEnsemble(
    structure=mc.structure, calculator=calc,
    temperature=temperature_min,
    forward=True,
    ensemble_data_write_interval=1,
    boltzmann_constant=1,
    n_steps=n_integration_steps)
mc.run()
data_container = mc.data_container

(_, free_energy_integration_forward) = \
        get_free_energy_thermodynamic_integration(data_container, cs,
                                                  forward=True,
                                                  max_temperature=temperature_max_plot_limit,
                                                  boltzmann_constant=1)

# Prepare for the TI by running the canonical ensemble at a low temperature.
mc = CanonicalEnsemble(
        structure=mc.structure,
        calculator=calc,
        temperature=temperature_min,
        boltzmann_constant=1,
        trajectory_write_interval=None,
        ensemble_data_write_interval=200)
mc.run(n_equilibration_steps)

mc = ThermodynamicIntegrationEnsemble(
    structure=mc.structure, calculator=calc,
    temperature=temperature_min,
    forward=False,
    ensemble_data_write_interval=1,
    boltzmann_constant=1,
    n_steps=n_integration_steps)
mc.run()
data_container = mc.data_container

(temperatures_integration, free_energy_integration_backward) = \
        get_free_energy_thermodynamic_integration(data_container, cs,
                                                  forward=False,
                                                  max_temperature=temperature_max_plot_limit,
                                                  boltzmann_constant=1)

# The best approximation for TI is now given by the average of the two runs.
free_energy_integration_average = 0.5 * (free_energy_integration_forward +
                                         free_energy_integration_backward)

##########################################################################################
# Temperature integration
##########################################################################################
# Prepare for the temperature integration by running the canonical ensemble at a high temperature.
mc = CanonicalEnsemble(
        structure=start_configuration,
        calculator=calc,
        temperature=temperature_max,
        boltzmann_constant=1,
        trajectory_write_interval=None,
        ensemble_data_write_interval=200)
mc.run(n_equilibration_steps)

mc = CanonicalAnnealing(
        structure=mc.structure,
        calculator=calc,
        T_start=temperature_max,
        T_stop=temperature_min,
        cooling_function='linear',
        n_steps=n_integration_steps,
        boltzmann_constant=1,
        trajectory_write_interval=None,
        ensemble_data_write_interval=1)
mc.run()
data_container = mc.data_container


data_container = mc.data_container
(temperatures_temperature, free_energy_temperature_forward) = \
    get_free_energy_temperature_integration(data_container,
                                            cs,
                                            forward=True,
                                            temperature_reference=temperature_max,
                                            max_temperature=temperature_max_plot_limit,
                                            boltzmann_constant=1)

# Prepare for the temperature integration by running the canonical ensemble at a low temperature.
mc = CanonicalEnsemble(
        structure=mc.structure,
        calculator=calc,
        temperature=temperature_min,
        boltzmann_constant=1,
        trajectory_write_interval=None,
        ensemble_data_write_interval=200)
mc.run(n_equilibration_steps)

mc = CanonicalAnnealing(
        structure=mc.structure,
        calculator=calc,
        T_start=temperature_min,
        T_stop=temperature_max,
        cooling_function='linear',
        n_steps=n_integration_steps,
        boltzmann_constant=1,
        trajectory_write_interval=None,
        ensemble_data_write_interval=1)
mc.run()
data_container = mc.data_container

(temperatures_temperature, free_energy_temperature_backward) = \
    get_free_energy_temperature_integration(data_container,
                                            cs,
                                            forward=False,
                                            temperature_reference=temperature_max,
                                            max_temperature=temperature_max_plot_limit,
                                            boltzmann_constant=1)

# The best approximation for TEI is now given by the average of the two runs.
free_energy_temperature_average = 0.5 * (free_energy_temperature_forward +
                                         free_energy_temperature_backward)

##########################################################################################
# Analytical solution
##########################################################################################
# Get all possible occupations
atoms = get_all_structures(supercell)

# Gets the free energy from the enumerated structures
temperatures = np.linspace(temperature_min, temperature_max_plot_limit, 1000)
free_energy_enumeration = get_free_energy_enumeration(temperatures,
                                                      atoms,
                                                      calc,
                                                      boltzmann_constant=1)

##########################################################################################
# Plotting
##########################################################################################
natoms = len(supercell)
fig, ax = plt.subplots()
ax.plot(temperatures_integration, free_energy_integration_forward / natoms, color='tab:red',
        label='forward', ls='-', lw=2)
ax.plot(temperatures_integration, free_energy_integration_backward / natoms, color='tab:blue',
        label='backward', ls='-', lw=2)
ax.plot(temperatures_integration, free_energy_integration_average / natoms, color='tab:grey',
        label='average', ls='-', lw=2)
ax.set_xlim((temperature_min, temperature_max_plot_limit))
ax.set_ylim((-45 / natoms, -30 / natoms))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('free energy / atoms (eV)')
ax.legend()
fig.tight_layout()
fig.savefig('thermodynamic_integration_free_energy_forward_backward.svg')

fig, ax = plt.subplots()
ax.plot(temperatures_temperature, free_energy_temperature_forward / natoms, color='tab:red',
        label='forward', ls='-', lw=2)
ax.plot(temperatures_temperature, free_energy_temperature_backward / natoms, color='tab:blue',
        label='backward', ls='-', lw=2)
ax.plot(temperatures_temperature, free_energy_temperature_average / natoms, color='tab:grey',
        label='average', ls='-', lw=2)
ax.set_xlim((temperature_min, temperature_max_plot_limit))
ax.set_ylim((-45 / natoms, -30 / natoms))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('free energy / atoms (eV)')
ax.legend()
fig.tight_layout()
fig.savefig('temperature_integration_free_energy_forward_backward.svg')

fig, ax = plt.subplots()
ax.axhline(-2, label='Ground state energy', color='black', ls='--')
ax.plot(temperatures, free_energy_enumeration / natoms, color='tab:red',
        label='enumeration', lw=2)
ax.plot(temperatures_integration, free_energy_integration_average / natoms, color='tab:blue',
        label='thermodynamic integration', ls='-', lw=2)
ax.plot(temperatures_temperature, free_energy_temperature_average / natoms, color='tab:grey',
        label='temperature integration', ls='-', lw=2)
ax.set_xlim((temperature_min, temperature_max_plot_limit))
ax.set_ylim((-45 / natoms, -30 / natoms))
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('free energy / atoms (eV)')
ax.legend()
fig.tight_layout()
fig.savefig('free_energy_canonical_ensemble.svg')
