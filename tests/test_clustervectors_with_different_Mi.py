"""
This scripts checks the computation of clustervectors with non constant Mi

"""

import numpy as np

import icetdev
from ase.build import bulk
from ase.db import connect
from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms

print(__doc__)


def test_mi_int_and_list(atoms, subelements, cutoffs):
    """Test that int Mi produces the same results as a list of Mi"""

    prim_size = len(
        icetdev.permutationMap.__get_primitive_structure(atoms.copy()))
    Mi_int = 2
    Mi_list = [2] * prim_size

    clusterspace_int = create_clusterspace(
        subelements, cutoffs, atoms=atoms, Mi=Mi_int)

    clusterspace_list = create_clusterspace(
        subelements, cutoffs, atoms=atoms, Mi=Mi_list)

    atoms_prim = clusterspace_int.get_primitive_structure().to_atoms()

    # create and populate a supercell and get cv
    repeat = [1] * 3
    for i, pbc in enumerate(atoms_prim.pbc):
        if pbc:
            repeat[i]=3
    conf = atoms_prim.repeat(repeat)
    for at in conf:
        at.symbol = np.random.choice(subelements)

    conf = structure_from_atoms(conf)

    cv_int = clusterspace_int.get_clustervector(conf)
    cv_list = clusterspace_list.get_clustervector(conf)
    assert (np.array(cv_int) == np.array(cv_list)).all, "Error: clustervectors were not equal for an equivalent Mi with: type(Mi) == int and type(Mi)==list"

db = connect("structures_for_testing.db")

subelements = ["H", "He", "Pb"]
for row in db.select():
    atoms_row = row.toatoms()
    atoms_tag = row.tag
    cutoffs = [1.4] * 3
    if len(atoms_row) == 0:
        continue
    if atoms_row.get_pbc().all() == True:
        atoms_row.wrap()
    if True:  # atoms_row.get_pbc().all() == True:
        print("Testing finding cv for structure: {} with cutoffs {}".format(
            atoms_tag, cutoffs))
        atoms_row.set_chemical_symbols(len(atoms_row) * [atoms_row[0].symbol])

        test_mi_int_and_list(atoms_row, subelements, cutoffs)
