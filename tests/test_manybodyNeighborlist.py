from icetdev import *
from icetdev.structure import *
from icetdev.manybodyNeighborlist import *
import numpy.random as random
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk
from tests import manybodyNeighborlistTester


neighbor_cutoff = 6.1

# set ut atoms and icet structure
a = bulk('Ti', "bcc", a=3.321).repeat(2)
a.set_pbc((True, True, True))
structure = structure_from_atoms(a)
# set up neighborlist for input to manybody neighborlist
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
neighbors = nl.get_neighbors(0)

# set up manybody neighborlist
mbnl = ManybodyNeighborlist()
# this is intersect between neighbors of atom 0 and atom 1
intersect = mbnl.calc_intersection(nl.get_neighbors(0), nl.get_neighbors(1))

# test intersect by doing a naive intersect
nbrs1 = nl.get_neighbors(0)
nbrs2 = nl.get_neighbors(1)
naive_intersect = []
for n1 in nbrs1:
    for n2 in nbrs2:
        if n1[0] == n2[0] and (n1[1] == n2[1]).all():
            naive_intersect.append(n1)


# assert that all the intersects are equal
for n1, n2 in zip(intersect, naive_intersect):
    assert n1[0] == n2[0] and (n1[1] == n2[1]).all()


# test actual mbnl
order = 3
bothways = True
index1 = 0
index2 = len(a) - 3
nbrs1 = mbnl.build(nl, index1, order, True)
nbrs2 = mbnl.build(nl, index2, order, True)
#print(len(nbrs1), len(nbrs2)) #debug
assert len(nbrs1) == len(
    nbrs2), "bothways = True should give same number of neigbhors independent on what index you look at"


# get manybodyNeighbors to third order
mbnl_T = manybodyNeighborlistTester.manybodyNeighborlistTester()
ase_nl = NeighborList(len(a) * [neighbor_cutoff / 2.0], skin=1e-8,
                      bothways=True, self_interaction=False)
ase_nl.update(a)


# nbrs = mbnl_T.build(ase_nl, index, order, bothways=bothways)
maxorder = 4
bothways = True
index = 0
for i in range(len(a)):
    for j in range(maxorder):
        index = i
        order = j
        nbrs_tester = mbnl_T.build(ase_nl, index, order, bothways)
        nbrs_cpp = mbnl.build(nl, index, order, bothways)
        assert len(nbrs_cpp) == len(nbrs_cpp)


bothways = False
for i in range(len(a)):
    for j in range(maxorder):
        index = j
        order = i
        nbrs_tester = mbnl_T.build(ase_nl, index, order, bothways)
        nbrs_cpp = mbnl.build(nl, index, order, bothways)
        assert len(nbrs_cpp) == len(nbrs_cpp)


# debug
def printNeighbor(nbr, onlyIndice=False):
    if onlyIndice:
        print(nbr[0], end=" ")
    else:
        print(nbr[0], nbr[1], end=" ")
