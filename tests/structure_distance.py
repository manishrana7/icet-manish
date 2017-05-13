from example import Structure
from ase import Atoms
import numpy as np
from ase.build import bulk

from icetdev import *
from icetdev.structure import *
from ase.neighborlist import NeighborList

#Distance tolerance for comparing distances 
DISTTOL = 1e-8

# Setup a chain of H,O,C
# H-O Dist = 2
# O-C Dist = 3
# C-H Dist = 5 with mic=False
# C-H Dist = 4 with mic=True
a = Atoms('HOC', positions=[(1, 1, 1), (3, 1, 1), (6, 1, 1)])
a.set_cell((9, 2, 2))
a.set_pbc((True, False, False))

ex_structure = structure_from_atoms(a)


# Calculate indiviually with mic=False
assert a.get_distance(0, 1, mic=False) == 2
assert a.get_distance(1, 2, mic=False) == 3
assert a.get_distance(0, 2, mic=False) == 5


# Test distance calculator with offsets
a = bulk("Al").repeat(2)
ex_structure = structure_from_atoms(a)
nl = NeighborList(len(a)*[10])
nl.update(a)
for ind in range(len(a)):
    indices, offsets = nl.get_neighbors(ind)
    for i, offset in zip(indices, offsets):
        dist = np.linalg.norm(a.positions[ind] - (a.positions[i] + np.dot(offset, a.get_cell())))
        dist_struct = ex_structure.get_distance2(ind,[0,0,0],i, offset)
        assert (dist - dist_struct) < DISTTOL