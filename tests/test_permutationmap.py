from icetdev.permutationMap import PermutationMap

from ase import Atoms
from ase.build import bulk

import spglib as spglib


#ASE atoms
atoms = bulk("Al","fcc",a=2.0).repeat([3,2,1])

symmetry = spglib.get_symmetry(atoms)
translations = symmetry['translations']
rotations = symmetry['rotations']
permutation_map = PermutationMap(translations, rotations)

pos = atoms.get_scaled_positions()
#print(pos)

assert len(rotations) == len(translations)
print(rotations)
print(translations)
print("pyrot size ", len(rotations))
permutation_map.build(pos)
perm_pos = permutation_map.get_permutated_positions()
ind_pos  = permutation_map.get_indiced_positions()
print(len(perm_pos))
for pp in perm_pos:
    print(pp)

for ind in list(map(list, zip(*ind_pos))):
    print(ind)

