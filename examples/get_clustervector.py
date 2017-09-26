from icetdev import orbitList
from icetdev.structure import structure_from_atoms, Structure
from icetdev.manybodyNeighborlist import *
from ase import Atoms
from ase.build import bulk
from icetdev.neighborlist import get_neighborlists
from icetdev.orbitList import create_orbit_list
from ase.spacegroup import crystal
import time
from icetdev.clusterCounts import *
from icetdev import ClusterSpace
import random
# from clib.cluster_space import ClusterSpace


def setup_test_orbitlist(atoms, cutoffs):
    structure = structure_from_atoms(atoms)
    mbnl = ManybodyNeighborlist()
    neighborlists = get_neighborlists(atoms=atoms, cutoffs=cutoffs)
    mbnl.build(neighborlists, 0, True)
    return structure, mbnl, neighborlists


latnbr = LatticeNeighbor(0, [0., 0., 0.])


atoms = bulk("Al", "bcc", a=1)
# atoms = crystal("Al", spacegroup=223)
from ase.io import read
# atoms = read(
#    "../clathrate-cluster-expansions/cluster-expansions/reference-clathrate-structure.json")
from ase.visualize import view
# view(atoms)
# exit(1)
cutoffs = [5,4,2.5,2.3] 

# cs = ClusterSpace(atoms=atoms, cutoffs=cutoffs)
# print(cs)
# len(cs)
structure, mbnl, neighborlists = setup_test_orbitlist(atoms, cutoffs)



############################################
#                                          #
#    orbitlist from permutation map        #
#                                          #
############################################

orbitlist = create_orbit_list(structure, cutoffs, verbosity=4)

# orbitlist.sort()

for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, orbitlist.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(orbitlist.size()))


#################################################
#                                               #
#             get clustervector                 #
#                                               #
#################################################

atoms_super = atoms.repeat(4)

for atom in atoms_super:
    if random.random() < 0.5:
        atom.symbol = 'He'

structure_atoms = structure_from_atoms(atoms_super)        

clusterspace = ClusterSpace(2,orbitlist)

t1=time.time()
cv1 = clusterspace.get_clustervector(structure_atoms)
t2=time.time()

print(cv1)

print("Time to get cv for supercell with {0} atoms and clusterspace of {1} clusters : {2:.4} s".format(len(atoms_super), orbitlist.size(), t2-t1))