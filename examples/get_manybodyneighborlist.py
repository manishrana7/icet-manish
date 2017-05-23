from icetdev import *
from icetdev.structure import *
from icetdev.manybodyNeighborlist import *
import numpy.random as random
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk


neighbor_cutoff = 5.0

# set ut atoms and icet structure
a = bulk("Al").repeat(2)
structure = structure_from_atoms(a)

# set up neighborlist for input to manybody neighborlist
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
neighbors = nl.get_neighbors(0)

# set up manybody neighborlist
mbnl = ManybodyNeighborlist()

#build mbnl for index "index" and order "order"
index = 0
order = 3
bothways = False
mbnl_i = mbnl.build(nl, index, order, bothways)
# get manybodyNeighbors to third order

def arrayCompare(arr1, arr2):
    assert len(arr1) == len(arr2)
    for i in range(3):
        if arr1[i] < arr2[i]:
            return True
        if arr1[i] > arr2[i]:
            return False
    return False

def nbrCompare(nbr1, nbr2):
    if nbr1[0] < nbr2[0]:
         return True
    if nbr1[0] > nbr2[0]:
         return False
    return arrayCompare(nbr1[1],nbr2[1])     
def naiveManybodyThirdOrder(nl, index, bothways=True):
    nbr_0 = nl.get_neighbors(index)
    nbr_index = (index, [0,0,0])
    nbrs = []
    for j in nbr_0:
        if not bothways and nbrCompare(j, nbr_index):
            continue
        nbr_j = nl.get_neighbors(j[0])
        for k in nbr_j:            
            if (nl.is_neighbor(index, k[0], k[1] + j[1])):
                neighbor_k = (k[0], k[1] + j[1])
                manybodyNbr = [index, [0., 0., 0.], j, neighbor_k]
                if not bothways:
                    if nbrCompare(j, neighbor_k):                    
                        nbrs.append(manybodyNbr)
                else:
                    nbrs.append(manybodyNbr)            
    return nbrs


def getDistances(nbrs, structure):
    d = []
    for i in nbrs:
        for j in nbrs:
            if nbrCompare(i,j) or i[0] == j[0] and (i[1] == j[1]).all():
                continue
            d.append( structure.get_distance2(i[0], i[1] ,j[0],j[1]))
    return d
# debug
def printNeighbor(nbr, onlyIndice=False):
    if onlyIndice:
        print(nbr[0], end=" ")
    else:
        print(nbr[0], nbr[1], end=" ")


mbnl_3 = naiveManybodyThirdOrder(nl, index, bothways)
for ss in mbnl_3:
    print(ss[0], ss[1], end=" ")
    printNeighbor(ss[2])
    printNeighbor(ss[3])
    print("")


for nbr in nl.get_neighbors(0):
    print(0, [0, 0, 0], end=" ")
    printNeighbor(nbr)
    print("")
for nbr in mbnl_i:
    orig, manyInd = nbr
    for mind in manyInd:
        nbr = []
        for s in orig:
            nbr.append(s)
            printNeighbor(s)
        printNeighbor(mind)
        nbr.append(mind)
        print(getDistances(nbr,structure))
        print("")

count = 0
for nbr in mbnl_i:
    orig, manyInd = nbr
    print("ec", len(orig), len(manyInd))
    count += len(manyInd)
print("counts", len(mbnl_3), count)

# Test if there are duplicates in mbnl
