from icetdev import orbitList
from icetdev.structure import structure_from_atoms, Structure
from icetdev.manybodyNeighborlist import *
from ase import Atoms
from ase.build import bulk
from icetdev.neighborlist import get_neighborlists
from icetdev.orbitList import create_orbit_list
from ase.spacegroup import crystal
import time
# from clib.cluster_space import ClusterSpace


def setup_test_orbitlist(atoms, cutoffs):
    structure = structure_from_atoms(atoms)
    mbnl = ManybodyNeighborlist()
    neighborlists = get_neighborlists(atoms=atoms, cutoffs=cutoffs)
    mbnl.build(neighborlists, 0, True)
    return structure, mbnl, neighborlists


atoms = bulk("Al", "bcc", a=1)
# atoms = crystal("Al", spacegroup=223)
from ase.io import read
atoms = read("../clathrate-cluster-expansions/cluster-expansions/reference-clathrate-structure.json")
from ase.visualize import view
# view(atoms)
# exit(1)
cutoffs = [5, 5.51]

# cs = ClusterSpace(atoms=atoms, cutoffs=cutoffs)
# print(cs)
# len(cs)
structure, mbnl, neighborlists = setup_test_orbitlist(atoms, cutoffs)

#ol = orbitList.OrbitList(mbnl, structure)
ol = orbitList.OrbitList(neighborlists, structure)

ol.sort()

#ol.print(verbosity = 3)
for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, ol.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(ol.size()))

# for orbit in range(ol.size()):
#     (ol.get_orbit(i).get_representative_cluster().print())
#     print("number of equivalent sites: ",len(ol.get_orbit(i).get_equivalent_sites()))

############################################
#                                          #
#   Test orbitlist from permutation map    #
#                                          #
############################################

# for nbr in neighborlists[0].get_neighbors(0):
#     print(nbr)
# print("-------")    
# for nbr in neighborlists[0].get_neighbors(1):
#     print(nbr)    
orbitlist = create_orbit_list(structure, cutoffs, verbosity=4)

#orbitlist.sort()

#ol.print(verbosity = 3)
for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, orbitlist.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(orbitlist.size()))

for i in range(orbitlist.size()):
     (orbitlist.get_orbit(i).get_representative_cluster().print())
     print("number of equivalent sites: ",len(orbitlist.get_orbit(i).get_equivalent_sites()))



#################################################
#                                               #
#        Get orbitlist for a supercell          #
#                                               #
#################################################


N = 3
atoms = atoms.repeat(N)

structure_repeat = structure_from_atoms(atoms)
t1 = time.time()
supercell_orbitlist = orbitlist.get_supercell_orbitlist(structure_repeat)
t2 = time.time()

print("Time to get supercell with x atoms {0}: {1} s".format(structure_repeat.size(), t2-t1))

print("size of repeated orbitlist {0}".format(supercell_orbitlist.size()))

for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, supercell_orbitlist.get_number_of_NClusters(i)))


for i in range(supercell_orbitlist.size()):
     (supercell_orbitlist.get_orbit(i).get_representative_cluster().print())
     print("number of equivalent sites: ",len(supercell_orbitlist.get_orbit(i).get_equivalent_sites()))


# import numpy as np
# for i in range(orbitlist.size()):
#     if len(orbitlist.get_orbit(i).get_representative_cluster().get_distances())==3:        
#         print(orbitlist.get_orbit(i).get_representative_cluster().get_distances())


#         eq_sites = orbitlist.get_orbit(i).get_equivalent_sites()
#         sites = eq_sites[0]
#         for s in sites:
#             print(np.dot(s.unitcellOffset, atoms.cell))
#         print(" end ")
#         print(" ")
# for i in range(orbitlist.size()):
#     if len(orbitlist.get_orbit(i).get_representative_cluster().get_distances())==3:        
#         dists = orbitlist.get_orbit(i).get_rep
# resentative_cluster().get_distances()
#         if np.linalg.norm((np.array(dists) -[0.86, 1.4,1.65])) < 0.1:
#             eq_sites = orbitlist.get_orbit(i).get_equivalent_sites()
#             print("len of eq sites: {}".format(len(eq_sites)) )
#             for sites in eq_sites:
#                 for s in sites:                    
#                     print(np.dot(s.unitcellOffset, atoms.cell), end= ' ')
#                 print( " ")                    
#             print(" end ")
#             print(" ")
