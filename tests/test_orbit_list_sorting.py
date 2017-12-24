from icetdev import Structure
from icetdev.orbit_list import create_orbit_list
from ase.build import bulk

'''
This test will construct two orbit lists each from the same atom structure but
one will have shuffled atoms. The test will assert that the orbits in the
orbit lists are the same from both orbit lists.

Todo
----
given that create_orbit_list get the primitive cell is this really a good test?
'''

atoms1 = bulk('Al', 'bcc', a=1).repeat(3)
cutoffs = [3, 2.5, 2]
structure = Structure.from_atoms(atoms1)


orbit_list_1 = create_orbit_list(structure, cutoffs, verbosity=0)

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

structure2 = Structure.from_atoms(atoms2)

# create orbit list from shuffled structure
orbit_list_2 = create_orbit_list(structure2, cutoffs, verbosity=4)


assert len(orbit_list_1) == len(orbit_list_2), \
    'The sizes of the different orbit lists should be equal'

orbit_list_1.sort()
orbit_list_2.sort()

for i in range(len(orbit_list_1)):
    assert orbit_list_1.get_orbit(i).get_representative_cluster() \
        == orbit_list_2.get_orbit(i).get_representative_cluster(), \
        'Representative clusters should be the same in the different orbits'