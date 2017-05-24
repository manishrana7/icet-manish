from tests import manybodyNeighborlistTester
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.build import bulk


mbnl_T = manybodyNeighborlistTester.manybodyNeighborlistTester()

atoms  = bulk("Al").repeat(2)

neighbor_cutoff = 5
ase_nl = NeighborList(len(atoms)*[neighbor_cutoff/2.0],skin=1e-8,
                    bothways=True,self_interaction=False)
ase_nl.update(atoms)


index = 1
order = 3
bothways = False


nbrs = mbnl_T.build(ase_nl, index, order, bothways=bothways) 

count = 0
for j in nbrs:    
      count += len(j[1])
      for intersect in j[1]:
            print(j[0], intersect)
print("count = {}".format(count))

