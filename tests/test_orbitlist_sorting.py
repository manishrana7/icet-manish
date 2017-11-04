from icetdev import orbitList, Structure
from icetdev.orbitList import create_orbit_list
from icetdev.structure import structure_from_atoms
from ase.build import bulk


# This test will construct two orbitlists each from the same atom structure but one will have shuffled atoms.
# The test will assert that the orbits in the orbitlists are the same from both orbitlists
# TODO: given that create_orbit_list get the primitive cell is this really
# a good test?

atoms1 = bulk(
    "Al", "bcc", a=1).repeat(3)
cutoffs = [3, 2.5, 2]
structure = structure_from_atoms(atoms1)


orbitlist_1 = create_orbit_list(structure, cutoffs, verbosity=0)

atoms2 = atoms1.copy()


# shuffle atoms2 to be different
for i in range(len(atoms2)):
    for j in range(len(atoms2)):
        if j < i:
            continue
        else:
            atom_temporary_hold = atoms2[i].position.copy()
            atoms2[i].position = atoms2[j].position.copy()
            atoms2[j].position = atom_temporary_hold


# just to make sure repeat the shuffled structure
atoms2 = atoms2.repeat(2)

structure2 = structure_from_atoms(atoms2)

#create orbitlist from shuffled structure
orbitlist_2 = create_orbit_list(structure2, cutoffs, verbosity=0)


assert len(orbitlist_1) == len(orbitlist_2), \
    "The sizes of the different orbitlists should be equal"

orbitlist_1.sort()
orbitlist_2.sort()

for i in range(len(orbitlist_1)):
    assert orbitlist_1.get_orbit(i).get_representative_cluster() == orbitlist_2.get_orbit(
        i).get_representative_cluster(), "Representative clusters should be the same in the different orbits"
