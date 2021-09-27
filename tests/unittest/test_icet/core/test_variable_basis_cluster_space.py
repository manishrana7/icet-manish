import pytest
from ase import Atom
from ase.build import bulk
from icet import VariableBasisClusterSpace


@pytest.fixture
def system(request):
    model = request.param

    # Create primitive structure and cluster space
    if model == 'binary_fcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='fcc', a=alat)
        cutoffs = [7, 6, 4]
        nparameters_per_orbit = [3] * 10 + [2] * 12
    elif model == 'ternary_fcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge', 'Ga']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='fcc', a=alat)
        cutoffs = [5, 4]
        nparameters_per_orbit = [0]
    elif model == 'binary_bcc':
        alat = 4.0
        chemical_symbols = [['Al', 'Ge']]
        prim = bulk(chemical_symbols[0][0], crystalstructure='bcc', a=alat)
        cutoffs = [10, 10, 4]
        nparameters_per_orbit = [4] * 7
    elif model == 'binary_hcp':
        alat, clat = 3.4, 5.1
        chemical_symbols = [['Ag', 'Pd'], ['Ag', 'Pd']]
        prim = bulk(chemical_symbols[0][0], a=alat, c=clat, crystalstructure='hcp')
        cutoffs = [5, 4]
        nparameters_per_orbit = [2] * 7
    elif model == 'ternary_hcp':
        alat, clat = 3.4, 5.1
        chemical_symbols = [['Ag', 'Pd', 'Cu'], ['Ag', 'Pd', 'Cu']]
        prim = bulk(chemical_symbols[0][0], a=alat, c=clat, crystalstructure='hcp')
        cutoffs = [5, 4]
        nparameters_per_orbit = [0]
    elif model == 'sublattices_fcc':
        alat = 4.0
        chemical_symbols = [['Ag', 'Pd'], ['H', 'X']]
        prim = bulk(chemical_symbols[0][0], a=alat, crystalstructure='fcc')
        prim.append(Atom('H', (alat / 2, alat / 2, alat / 2)))
        cutoffs = [5, 4]
        nparameters_per_orbit = [0]
    elif model == 'inactivesublattice_fcc':
        alat = 4.0
        chemical_symbols = [['Pd', 'Ag'], ['W']]
        prim = bulk(chemical_symbols[0][0], a=alat, crystalstructure='fcc')
        prim.append(Atom('W', (alat / 2, alat / 2, alat / 2)))
        cutoffs = [5, 4]
        nparameters_per_orbit = [2] * 6
    elif model == 'inactivesublatticebutsamespecies_fcc':
        alat = 4.0
        chemical_symbols = [['Ag', 'Pd'], ['Ag']]
        prim = bulk(chemical_symbols[0][0], a=alat, crystalstructure='fcc')
        prim.append(Atom('Ag', (alat / 2, alat / 2, alat / 2)))
        cutoffs = [5, 4]
        nparameters_per_orbit = [2] * 6
    else:
        raise Exception(f'Unknown model ({model})')

    return prim, cutoffs, nparameters_per_orbit, chemical_symbols


def get_complementary_chemical_symbol(symbol, chemical_symbols):
    for site_symbols in chemical_symbols:
        if symbol == site_symbols[0]:
            return site_symbols[-1]
    raise Exception('Failed to find complementary chemcial symbol')


def get_supercell(primitive_structure, chemical_symbols, repeat, pattern):
    structure = primitive_structure.repeat(repeat)
    anti_structure = primitive_structure.repeat(repeat)
    for atom in anti_structure:
        atom.symbol = get_complementary_chemical_symbol(atom.symbol, chemical_symbols)
    anti_structure = anti_structure.repeat(repeat)

    # Occupy supercell according to some pattern
    if pattern == 'homogeneous':
        pass
    elif pattern == 'pseudorandom':
        for i in [2, 3, 4, 7, 11, 14, 15, 16, 17]:
            if i >= len(structure):
                break
            current_symbol = structure[i].symbol
            structure[i].symbol = anti_structure[i].symbol
            anti_structure[i].symbol = current_symbol
        if 'ternary' in pattern:
            structure[0].symbol = chemical_symbols[0][1]
    elif pattern == 'ordered':
        for i in range(len(structure)):
            if i % 2 == 1:
                continue
            current_symbol = structure[i].symbol
            if 'ternary' in pattern and i % 4 == 0:
                structure[i].symbol = chemical_symbols[0][1]
            else:
                structure[i].symbol = anti_structure[i].symbol
            anti_structure[i].symbol = current_symbol
    elif pattern == 'segregated':
        for i in range(len(structure) // 2):
            current_symbol = structure[i].symbol
            if 'ternary' in pattern:
                if i % 2 == 0:
                    structure[i].symbol = chemical_symbols[0][1]
                else:
                    structure[i].symbol = anti_structure[i].symbol
            else:
                structure[i].symbol = anti_structure[i].symbol
            anti_structure[i].symbol = current_symbol

    else:
        raise Exception(f'Unknown supercell ({structure})')

    return structure


@pytest.mark.parametrize('system, expected_active_symbols, expected_concentration_symbol', [
    ('binary_fcc', ('Al', 'Ge'), 'Ge'),
    ('binary_hcp', ('Ag', 'Pd'), 'Ag'),
    ('inactivesublattice_fcc', ('Ag', 'Pd'), 'Ag'),
    ('inactivesublatticebutsamespecies_fcc', ('Ag', 'Pd'), 'Ag')],
    indirect=['system'])
def test_initialization(system, expected_active_symbols, expected_concentration_symbol):
    prim, cutoffs, nparameters_per_orbit, symbols = system
    vbcs = VariableBasisClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=symbols,
                                     nparameters_per_orbit=nparameters_per_orbit)
    assert isinstance(vbcs, VariableBasisClusterSpace)
    assert vbcs.active_chemical_symbols == expected_active_symbols
    assert vbcs.concentration_symbol == expected_concentration_symbol


@pytest.mark.parametrize('system, expectation', [
    ('ternary_fcc', pytest.raises(NotImplementedError)),
    ('ternary_hcp', pytest.raises(NotImplementedError)),
    ('sublattices_fcc', pytest.raises(NotImplementedError))],
    indirect=['system'])
def test_initialization_for_systems_not_handled(system, expectation):
    """
    Test that errors are raised when initializating with systems
    that are not handled.
    """
    prim, cutoffs, nparameters_per_orbit, symbols = system
    with expectation:
        _ = VariableBasisClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=symbols,
                                      nparameters_per_orbit=nparameters_per_orbit)


@pytest.mark.parametrize('system, nparameters_per_orbit', [
    ('binary_fcc', [1, 2, 3])],
    indirect=['system'])
def test_initialization_with_incommensurate_nparameters_per_orbit(system, nparameters_per_orbit):
    """
    Test that an error is raised when nparameters_per_orbit
    does not match the cluster space.
    """
    prim, cutoffs, _, symbols = system
    with pytest.raises(ValueError):
        _ = VariableBasisClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=symbols,
                                      nparameters_per_orbit=nparameters_per_orbit)


@pytest.mark.parametrize('system, repeat, pattern, expected_concentration', [
    ('binary_fcc', (3, 3, 3), 'pseudorandom', 9 / 27),
    ('binary_hcp', (3, 3, 3), 'pseudorandom', 45 / 54),
    ('inactivesublattice_fcc', (3, 3, 3), 'pseudorandom', 4 / 27),
    ('inactivesublatticebutsamespecies_fcc', (4, 3, 3), 'pseudorandom', 32 / 36)],
    indirect=['system'])
def test_get_concentration(system, repeat, pattern, expected_concentration):
    prim, cutoffs, nparameters_per_orbit, symbols = system
    supercell = get_supercell(prim, symbols, repeat, pattern)
    vbcs = VariableBasisClusterSpace(prim, cutoffs=cutoffs, chemical_symbols=symbols,
                                     nparameters_per_orbit=nparameters_per_orbit)
    retval = vbcs.get_concentration(supercell)
    assert abs(retval - expected_concentration) < 1e-6
