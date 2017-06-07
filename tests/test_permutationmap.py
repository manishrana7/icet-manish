from icetdev.permutationMap import PermutationMap

from ase import Atoms
from ase.build import bulk

import spglib as spglib


#ASE atoms
atoms = bulk("Al","fcc",a=2).repeat(3)

symmetry = spglib.get_symmetry(atoms)
print(symmetry['translations'])
print(symmetry['rotations'])
#permutation_map = PermutationMap(symmetry['translations'], symmetry['rotations'])
