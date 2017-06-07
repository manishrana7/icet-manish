from icetdev.permutationMap import PermutationMap

from ase import Atoms
from ase.build import bulk

import spglib as spglib


#ASE atoms
atoms = bulk("Al","fcc",a=2).repeat(3)

symmetry = spglib.get_symmetry(atoms)
translations = symmetry['translations']
rotations = symmetry['rotations']
permutation_map = PermutationMap(translations, rotations)

pos = atoms.get_scaled_positions()
#print(pos)

assert len(rotations) == len(translations)
permutation_map.build(pos)
