import pytest

from ase import Atom
from ase.build import bulk
from icet import ClusterExpansion, ClusterSpace
from mchammer.calculators.cluster_expansion_calculator import \
    ClusterExpansionCalculator
import numpy as np


def get_complementary_symbol(symbol, chemical_symbols):
    for site_symbols in chemical_symbols:
        if symbol == site_symbols[0]:
            return site_symbols[-1]
    raise Exception('Failed finding complementary symbol')


@pytest.fixture
def system(request):
    model, repeat, supercell = request.param

    # Create primitive structure and cluster space
    if model == 'binary_fcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='fcc', a=alat)
        cutoffs = [5, 5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'ternary_fcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge', 'Ga']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='fcc', a=alat)
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'binary_bcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='bcc', a=alat)
        cutoffs = [10, 10, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'ternary_bcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge', 'Ga']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='bcc', a=alat)
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'binary_hcp':
        alat, clat = 3.4, 5.1
        chemical_symbols = [['Ag', 'Pd'], ['Ag', 'Pd']]
        prim = bulk(chemical_symbols[0][0], a=alat, c=clat, crystalstructure='hcp')
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'ternary_hcp':
        alat, clat = 3.4, 5.1
        chemical_symbols = [['Ag', 'Pd', 'Cu'], ['Ag', 'Pd', 'Cu']]
        prim = bulk(chemical_symbols[0][0], a=alat, c=clat, crystalstructure='hcp')
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'sublattices_fcc':
        alat = 4.0
        chemical_symbols = [['Ag', 'Pd'], ['H', 'X']]
        prim = bulk(chemical_symbols[0][0], a=alat, crystalstructure='fcc')
        prim.append(Atom('H', (alat / 2, alat / 2, alat / 2)))
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'ternarysublattices_fcc':
        alat = 4.0
        chemical_symbols = [['Ag', 'Pd', 'Cu'], ['H', 'X']]
        prim = bulk(chemical_symbols[0][0], a=alat, crystalstructure='fcc')
        prim.append(Atom('H', (alat / 2, alat / 2, alat / 2)))
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    elif model == 'inactivesublattice_fcc':
        alat = 4.0
        chemical_symbols = [['Ag', 'Pd'], ['W']]
        prim = bulk(chemical_symbols[0][0], a=alat, crystalstructure='fcc')
        prim.append(Atom('W', (alat / 2, alat / 2, alat / 2)))
        cutoffs = [5, 4]
        cs = ClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=chemical_symbols)
    else:
        raise Exception(f'Unknown model ({model})')

    # Make a supercell as well as an "inverted" version of the supercell,
    # the latter with the purpose of more easily decide how atoms can be changed
    structure = prim.repeat((2, 2, 3))
    anti_structure = prim.copy()
    for atom in anti_structure:
        atom.symbol = get_complementary_symbol(atom.symbol, chemical_symbols)
    anti_structure = anti_structure.repeat((2, 2, 3))

    # Occupy supercell according to some pattern
    if supercell == 'homogeneous':
        pass
    elif supercell == 'pseudorandom':
        for i in [2, 3, 4, 7, 11, 14, 15, 16, 17]:
            if i >= len(structure):
                break
            current_symbol = structure[i].symbol
            structure[i].symbol = anti_structure[i].symbol
            anti_structure[i].symbol = current_symbol
        if 'ternary' in model:
            structure[0].symbol = chemical_symbols[0][1]
    elif supercell == 'ordered':
        for i in range(len(structure)):
            if i % 2 == 1:
                continue
            current_symbol = structure[i].symbol
            if 'ternary' in model and i % 4 == 0:
                structure[i].symbol = chemical_symbols[0][1]
            else:
                structure[i].symbol = anti_structure[i].symbol
            anti_structure[i].symbol = current_symbol
    elif supercell == 'segregated':
        for i in range(len(structure) // 2):
            current_symbol = structure[i].symbol
            if 'ternary' in model:
                if i % 2 == 0:
                    structure[i].symbol = chemical_symbols[0][1]
                else:
                    structure[i].symbol = anti_structure[i].symbol
            else:
                structure[i].symbol = anti_structure[i].symbol
            anti_structure[i].symbol = current_symbol

    else:
        raise Exception(f'Unknown supercell ({supercell})')

    # Define ECIs that are not all the same
    params = [(-1)**i * ((i + 1) / 10)**1.02 for i in range(len(cs))]
    ce = ClusterExpansion(cluster_space=cs, parameters=params)
    return ce, structure, anti_structure


# Make a list of parameters; possible combinations of systems and supercells
systems = []
systems_with_calculator_choice = []
for model in ['binary_fcc', 'ternary_fcc', 'binary_bcc', 'ternary_hcp',
              'sublattices_fcc', 'ternarysublattices_fcc', 'inactivesublattice_fcc']:
    for repeat in [(1, 1, 1), (2, 1, 1), (2, 2, 3)]:
        for supercell in ['homogeneous', 'pseudorandom', 'ordered', 'segregated']:
            if repeat in [(1, 1, 1), (2, 1, 1)] and supercell != 'ordered':
                continue
            elif 'ternary' in model and supercell == 'segregated':
                continue
            systems.append(((model, repeat, supercell)))
            systems_with_calculator_choice.append(((model, repeat, supercell), True))
            if model == 'binary_bcc':
                systems_with_calculator_choice.append(((model, repeat, supercell), False))


@pytest.mark.parametrize('system', systems, indirect=['system'])
def test_initialization(system):
    ce, structure, _ = system
    calc = ClusterExpansionCalculator(structure, ce, name='Test CE calc')
    assert isinstance(calc, ClusterExpansionCalculator)
    assert isinstance(calc.cluster_expansion, ClusterExpansion)
    assert calc.name == 'Test CE calc'
    assert abs(calc._property_scaling - len(structure)) < 1e-6
    assert calc.use_local_energy_calculator

    # Some alternative input parameters
    calc = ClusterExpansionCalculator(structure, ce, scaling=5.0, use_local_energy_calculator=False)
    assert isinstance(calc, ClusterExpansionCalculator)
    assert isinstance(calc.cluster_expansion, ClusterExpansion)
    assert calc.name == 'Cluster Expansion Calculator'
    assert abs(calc._property_scaling - 5.0) < 1e-6
    assert not calc.use_local_energy_calculator


@pytest.mark.parametrize('system', systems, indirect=['system'])
def test_get_cluster_vector(system):
    """Tests retrieval of full cluster vector from C++ side calculator against
    full cluster vector calculation from cluster space."""
    ce, structure, anti_structure = system
    calc = ClusterExpansionCalculator(structure, ce, name='Test CE calc')
    cv_calc = calc.cpp_calc.get_cluster_vector(structure.get_atomic_numbers())
    cv_cs = ce.get_cluster_space_copy().get_cluster_vector(structure)
    assert np.allclose(cv_calc, cv_cs)

    # Make sure it works after modifying the structure
    for i in range(2):
        structure[i].symbol = anti_structure[i].symbol
    cv_calc = calc.cpp_calc.get_cluster_vector(structure.get_atomic_numbers())
    cv_cs = ce.get_cluster_space_copy().get_cluster_vector(structure)
    assert np.allclose(cv_calc, cv_cs)


@pytest.mark.parametrize('system, use_local_energy_calculator',
                         systems_with_calculator_choice[:30],
                         indirect=['system'])
def test_change_calculation_flip(system, use_local_energy_calculator):
    """Tests differences when flipping."""
    ce, structure, anti_structure = system
    calc = ClusterExpansionCalculator(structure, ce, name='Test CE calc',
                                      use_local_energy_calculator=use_local_energy_calculator)
    for i in range(len(structure)):
        if structure[i].symbol == 'W':
            # inactive site
            continue
        sites = [i]
        change_local, change_global, change_ce = \
            get_energy_changes(calc, structure, anti_structure, sites)
        assert abs(change_local - change_global) < 1e-6
        assert abs(change_global - change_ce) < 1e-6


@pytest.mark.parametrize('system, use_local_energy_calculator',
                         systems_with_calculator_choice,
                         indirect=['system'])
def test_change_calculation_swap(system, use_local_energy_calculator):
    """Tests differences when swapping."""
    ce, structure, anti_structure = system
    calc = ClusterExpansionCalculator(structure, ce, name='Test CE calc',
                                      use_local_energy_calculator=use_local_energy_calculator)
    print(structure, anti_structure)
    for i in range(len(structure)):
        for j in range(3):
            if structure[i].symbol == 'W' or structure[j].symbol == 'W':
                # inactive site
                continue
            if j >= len(structure) or i == j:
                continue
            print(i, j, structure)
            sites = [i, j]
            change_local, change_global, change_ce = \
                get_energy_changes(calc, structure, anti_structure, sites)
            print(change_local, change_global)
            assert abs(change_local - change_global) < 1e-6
            assert abs(change_global - change_ce) < 1e-6


def get_energy_changes(calc, structure, anti_structure, sites):
    """
    Calculates change in property upon some change in a structure
    with three different methods.
    """
    change_local = calc.calculate_change(
        sites=sites,
        current_occupations=structure.get_atomic_numbers(),
        new_site_occupations=anti_structure.get_atomic_numbers()[sites])

    occupations_before = structure.get_atomic_numbers()
    occupations_after = structure.get_atomic_numbers()
    occupations_after[sites] = anti_structure.get_atomic_numbers()[sites]
    change_global = calc.calculate_total(occupations=occupations_after) \
        - calc.calculate_total(occupations=occupations_before)

    e_before_ce = calc.cluster_expansion.predict(structure)
    structure_copy = structure.copy()
    structure_copy.set_atomic_numbers(occupations_after)
    e_after_ce = calc.cluster_expansion.predict(structure_copy)
    change_ce = e_after_ce - e_before_ce
    change_ce *= len(structure)

    return change_local, change_global, change_ce
