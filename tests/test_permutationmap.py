from icetdev.permutationMap import PermutationMap

from ase import Atoms
from ase.build import bulk

import spglib as spglib


#ASE atoms
atoms = bulk("Al","fcc",a=2.0).repeat(1)

symmetry = spglib.get_symmetry(atoms)
translations = symmetry['translations']
rotations = symmetry['rotations']
permutation_map = PermutationMap(translations, rotations)

pos = atoms.get_scaled_positions()
#print(pos)

assert len(rotations) == len(translations)
print(rotations)
print("pyrot size ", len(rotations))
permutation_map.build(pos)

print( permutation_map.get_permutated_positions() )