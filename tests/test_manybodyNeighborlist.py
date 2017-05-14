from icetdev import *
from icetdev.structure import *
from icetdev.manybodyNeighborlist import *
import numpy.random as random
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk


neighbor_cutoff = 5

#set ut atoms and icet structure
a = bulk('Ti',"bcc",a=3.321).repeat(2)
a.set_pbc((True, True, True))
structure = structure_from_atoms(a)
#set up neighborlist for input to manybody neighborlist
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
neighbors = nl.get_neighbors(0)

#set up manybody neighborlist
mbnl = ManybodyNeighborlist()
#this is intersect between neighbors of atom 0 and atom 1
intersect = mbnl.calc_intersection(nl.get_neighbors(0), nl.get_neighbors(1))

#test intersect by doing a naive intersect
nbrs1 = nl.get_neighbors(0)
nbrs2 = nl.get_neighbors(1)
naive_intersect = []
for n1 in nbrs1:    
    for n2 in nbrs2:
        if n1[0] == n2[0] and (n1[1] == n2[1]).all():
            naive_intersect.append(n1)
"""
debug    
print("")
for n1 in naive_intersect:
    print(n1[0],n1[1][0],n1[1][1],n1[1][2])
"""
#assert that all the intersects are equal
for n1, n2 in zip(intersect, naive_intersect):
    assert n1[0] == n2[0] and (n1[1] == n2[1]).all()