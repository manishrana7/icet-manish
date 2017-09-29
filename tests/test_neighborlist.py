from icetdev import *
from icetdev.structure import *
import numpy.random as random
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk

neighbor_cutoff = 6

#ASE neighborlist
atoms = bulk("Al","fcc",a=2).repeat(3)
atoms = bulk('Ti', "bcc", a=3.321).repeat(2)

atoms.pbc = [True,False,True]
structure = structure_from_atoms(atoms)
#skin here should be thought of as DISTTOL in icet
ase_nl = NeighborList(len(atoms)*[neighbor_cutoff/2],skin=1e-8,
                    bothways=True,self_interaction=False)
ase_nl.update(atoms)
ase_indices, ase_offsets = ase_nl.get_neighbors(1)

# icet neighborlist
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
neighbors = nl.get_neighbors(1)
indices = []
offsets = []
for nbr in neighbors:
    indices.append(nbr.index)
    offsets.append(nbr.unitcellOffset)
    
"""
For debugging
print(len(indices), len(ase_indices))
for ind, offset in zip(indices, offsets):
    print(ind, offset, structure.get_distance2(0,[0,0,0],ind,offset))
print("====================================")
for ind, offset in zip(ase_indices, ase_offsets):
    print(ind, offset, structure.get_distance2(0,[0,0,0],ind,offset))

"""

assert len(indices) == len(ase_indices)
assert len(ase_offsets) == len(offsets)


#assert that each offset in icet neighborlist exists in ASE version
#assert that there are no duplicate offsets
#assert that the index for the corresponding offset are eq_indices 
for i,offset in zip(indices, offsets):
    assert offset in ase_offsets    
    eq_indices = [x for x, ase_offset in enumerate( ase_offsets) if ase_indices[x] == i and (ase_offset == offset).all()]
    if len(eq_indices) > 1:
        print(i, offset, eq_indices)
    assert len(eq_indices) == 1
    assert i == ase_indices[ eq_indices[0] ]

#test that the number of neighbors of  first and last is equal.
first_neighbors = nl.get_neighbors(0)
last_neighbors = nl.get_neighbors(len(atoms)-1)
assert len(first_neighbors) == len(last_neighbors)
