from icetdev import *
from icetdev.structure import *
from icetdev.manybodyNeighborlist import *
import numpy.random as random
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk


neighbor_cutoff = 5.1

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

bothways = True
mbnl_i = mbnl.build(nl,0,3,bothways)
#get manybodyNeighbors to third order
def naiveManybodyThirdOrder(nl, index, bothways=True):
    nbr_0 = nl.get_neighbors(index)
    nbrs = []
    for j in nbr_0:
        nbr_j = nl.get_neighbors(j[0])
        for k in nbr_j:

            if(not bothways and k[0] < j[0] and (k[1] == [0,0,0]).all() and (k[1] == j[1]).all()):
                continue
            

            if (nl.is_neighbor(index,k[0],k[1] + j[1])):
                neighbor_k = (k[0],k[1] + j[1])                
                manybodyNbr = [index,[0.,0.,0.],j,neighbor_k]
                nbrs.append(manybodyNbr)
    return nbrs                        


#debug
def printNeighbor(nbr,onlyIndice=False):
    if onlyIndice:
        print(nbr[0], end=" ")    
    else:        
        print(nbr[0], nbr[1], end=" ")

mbnl_3 = naiveManybodyThirdOrder(nl,0,bothways)
for ss in mbnl_3:
    print(ss[0],ss[1], end= " ")
    printNeighbor(ss[2])
    printNeighbor(ss[3])
    print("")



for nbr in nl.get_neighbors(0):
    print(0,[0,0,0], end= " ")
    printNeighbor(nbr)
    print("")
for nbr in mbnl_i:
    orig, manyInd = nbr
    for mind in manyInd:
        for s in orig:
            printNeighbor(s)                
        printNeighbor(mind)
        print("")    

count = 0
for nbr in mbnl_i:
    orig, manyInd = nbr    
    print("ec",len(orig),len(manyInd))
    count += len(manyInd)
print("counts",len(mbnl_3), count)

#Test if there are duplicates in mbnl

