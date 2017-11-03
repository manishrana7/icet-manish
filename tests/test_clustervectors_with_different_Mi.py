'''
This scripts checks the computation of cluster vectors with a variable number
of components on different sites.
'''

import numpy as np

import icetdev
from ase.db import connect
from icetdev.clusterspace import create_clusterspace, get_singlet_info
from icetdev.structure import structure_from_atoms

print(__doc__)


def test_mi_int_list_and_dict(atoms, subelements, cutoffs, allowed_sites):
    '''
    Test that int Mi, list Mi and dict Mi produces the same clustervectors.
    '''

    prim_size = len(
        icetdev.permutationMap.__get_primitive_structure(atoms.copy()))

    Mi_int = allowed_sites
    Mi_list = [allowed_sites] * prim_size
    Mi_dict = {}

    singlet_data = get_singlet_info(atoms.copy())
    for singlet in singlet_data:
        Mi_dict[singlet['orbit_index']] = allowed_sites

    clusterspace_int = create_clusterspace(atoms, cutoffs, subelements,
                                           Mi=Mi_int)

    clusterspace_list = create_clusterspace(atoms, cutoffs, subelements,
                                            Mi=Mi_list)

    clusterspace_dict = create_clusterspace(atoms, cutoffs, subelements,
                                            Mi=Mi_dict)

    atoms_prim = clusterspace_int.get_primitive_structure().to_atoms()

    # create and populate a supercell and get cv
    repeat = [1] * 3
    for i, pbc in enumerate(atoms_prim.pbc):
        if pbc:
            repeat[i] = 3
    conf = atoms_prim.repeat(repeat)
    for at in conf:
        at.symbol = np.random.choice(subelements)

    conf = structure_from_atoms(conf)

    cv_int = clusterspace_int.get_clustervector(conf)
    cv_list = clusterspace_list.get_clustervector(conf)
    cv_dict = clusterspace_dict.get_clustervector(conf)

    msg = 'Cluster vectors were not equal for an equivalent Mi (int vs list)'
    assert (np.array(cv_int) == np.array(cv_list)).all, msg

    msg = 'Cluster vectors were not equal for an equivalent Mi (int vs dict)'
    assert (np.array(cv_int) == np.array(cv_dict)).all, msg


db = connect('structures_for_testing.db')

subelements = ['H', 'He', 'Pb']
for row in db.select():
    atoms_row = row.toatoms()
    atoms_tag = row.tag
    cutoffs = [1.4] * 3
    if len(atoms_row) == 0:
        continue
    if atoms_row.get_pbc().all():
        atoms_row.wrap()
    if True:  # atoms_row.get_pbc().all() == True:
        print('Testing finding cv for structure: {} with cutoffs {}'.format(
            atoms_tag, cutoffs))
        atoms_row.set_chemical_symbols(len(atoms_row) * [atoms_row[0].symbol])
        for allowed_sites in range(2, 4):
            test_mi_int_list_and_dict(atoms_row, subelements,
                                      cutoffs, allowed_sites)
