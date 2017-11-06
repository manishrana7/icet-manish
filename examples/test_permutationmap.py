from icetdev.permutation_map import PermutationMap
from icetdev.tools.geometry import get_scaled_positions

from icetdev.structure import *
from icetdev.neighborlist import *
from ase import Atoms
from ase.build import bulk

import spglib as spglib

import numpy as np
# ASE atoms
atoms = bulk("Al", "fcc", a=2.0).repeat(1)

#get neighborlist
neighbor_cutoff = 2.0 #Ångström
structure = structure_from_atoms(atoms)
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
pos_neighbors = []
for latNbr in nl.get_neighbors(0):
    pos = structure.get_position(latNbr)
    pos_neighbors.append(pos)

frac_coordinates = get_scaled_positions(np.array(pos_neighbors), cell=atoms.cell, wrap=False, pbc = structure.pbc)

for fpos, pos in zip(frac_coordinates,pos_neighbors):
    print(fpos, pos, np.dot(fpos, atoms.cell), pos - np.dot(fpos, atoms.cell))

#print(structure.cell)
#exit(1)

symmetry = spglib.get_symmetry(atoms)
translations = symmetry['translations']
rotations = symmetry['rotations']
permutation_map = PermutationMap(translations, rotations)

#pos = atoms.get_scaled_positions()


assert len(rotations) == len(translations)

print("len of pos {}".format(len(frac_coordinates)))
permutation_map.build(frac_coordinates)
perm_pos = permutation_map.get_permutated_positions()
print("size of perm pos",len(perm_pos))
ind_pos, unique_pos = permutation_map.get_indiced_positions()
#exit(1)
print("Permutated fractional coordinates")
print("len of permutation pos",len(perm_pos))
for pp in perm_pos:
    unique_rows = np.vstack({tuple(row) for row in pp})
    for el in unique_rows:
        print(el, end=' ')
    print("")        
    

print("indices positions")

for i,pos in enumerate(ind_pos):
    print(i,len(set(pos)), pos)
#for ind in list(map(list, zip(*ind_pos))):
#    print(ind[0], set(ind))

print("index and unique position")
for index, dist in enumerate(unique_pos):
    print(index, dist)
