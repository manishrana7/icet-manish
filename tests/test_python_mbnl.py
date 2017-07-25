from tests import manybodyNeighborlistTester
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk


mbnl_T = manybodyNeighborlistTester.manybodyNeighborlistTester()

atoms  = bulk("Al").repeat(2)

neighbor_cutoff = 6.3

# set ut atoms and icet structure
atoms = bulk('Ti', "bcc", a=3.321).repeat(2)
atoms.set_pbc((True, True, True))


ase_nl = NeighborList(len(atoms)*[neighbor_cutoff/2.0],skin=1e-8,
                    bothways=True,self_interaction=False)
ase_nl.update(atoms)


index = 1
order = 3
bothways = False


nbrs = mbnl_T.build(ase_nl, index, order, bothways=bothways) 



# count = 0
# for j in nbrs:    
#       count += len(j[1])
#       for intersect in j[1]:
#             print(j[0], intersect)
# print("count = {}".format(count))
mbnl_T = manybodyNeighborlistTester.manybodyNeighborlistTester()

order = 3
bothways = True
index1 = 0
index2 = len(atoms) - 1

nbrs1 = mbnl_T.build(ase_nl, index1, order, True)
nbrs2 = mbnl_T.build(ase_nl, index2, order, True)
#print(len(nbrs1), len(nbrs2)) #debug
assert len(nbrs1) == len(
    nbrs2), "bothways = True should give same number of neigbhors independent on what index you look at. {} != {}".format(len(nbrs1), len(nbrs2))

print(len(nbrs1), len(nbrs2))