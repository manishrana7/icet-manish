import numpy as np
import numpy.random as random

from ase import Atoms
from ase.build import bulk
from ase.neighborlist import NeighborList
from icetdev import *
from icetdev.manybodyNeighborlist import *
from icetdev.structure import *
from tests import manybodyNeighborlistTester

# note that currently test at row 44 fails if cutoff is 6.1
neighbor_cutoff = 6.1

# set ut atoms and icet structure
a = bulk('Ti', "bcc", a=3.321).repeat(3)
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
        if n1.index == n2.index and (n1.unitcellOffset == n2.unitcellOffset).all():
            naive_intersect.append(n1)


# assert that all the intersects are equal
for n1, n2 in zip(intersect, naive_intersect):
    assert n1.index == n2.index and (
        n1.unitcellOffset == n2.unitcellOffset).all()


# test actual mbnl
order = 5
bothways = True
index1 = 0
index2 = len(a) - 1
nbrs1 = mbnl.build(order * [nl], index1, True)
nbrs2 = mbnl.build(order * [nl], index2, True)
# print(nbrs1, nbrs2) #debug

assert len(nbrs1) == len(
    nbrs2), "bothways = True should give same number of neigbhors independent on what index you look at. {} != {}".format(len(nbrs1), len(nbrs2))


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
    for j in range(2, maxorder):
        index = i
        order = j
        nbrs_tester = mbnl_T.build(
            (order - 1) * [ase_nl], index, bothways)
        nbrs_cpp = mbnl.build((order - 1) * [nl], index, bothways)
        assert len(nbrs_tester) == len(nbrs_cpp), "Number of manybody-neighbors at "\
            "index {0} with order {1} was not equal. {2} != {3}".format(
                i, j, len(nbrs_tester), len(nbrs_cpp))


# test that bothways = false also works
bothways = False
for i in range(len(a)):
    for j in range(1, maxorder):
        index = i
        order = j
        nbrs_tester = mbnl_T.build(order * [ase_nl], index, bothways)
        nbrs_cpp = mbnl.build(order * [nl], index, bothways)
        assert len(nbrs_tester) == len(nbrs_cpp), "python mbnl and cpp mbnl do not give same amount of "\
            "neighbors at index {0} with order {1} was not equal. {2} != {3}".format(
            i, j, len(nbrs_tester), len(nbrs_cpp))


# debug
def printNeighbor(nbr, onlyIndice=False):
    if onlyIndice:
        print(nbr.index, end=" ")
    else:
        print(nbr.index, nbr.unitcellOffset, end=" ")
