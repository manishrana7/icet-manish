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
cutoffs = [1.61] * 4

# cs = ClusterSpace(atoms=atoms, cutoffs=cutoffs)
# print(cs)
# len(cs)
structure, mbnl, neighborlists = setup_test_orbitlist(atoms, cutoffs)

# ol = orbitList.OrbitList(mbnl, structure)
orbitlist_no_symmetry = orbitList.OrbitList(neighborlists, structure)

orbitlist_no_symmetry.sort()

# ol.print(verbosity = 3)
for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, orbitlist_no_symmetry.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(orbitlist_no_symmetry.size()))

# for orbit in range(ol.size()):
#     (ol.get_orbit(i).get_representative_cluster().print())
# print("number of equivalent sites:
# ",len(ol.get_orbit(i).get_equivalent_sites()))

############################################
#                                          #
#   Test orbitlist from permutation map    #
#                                          #
############################################

orbitlist = create_orbit_list(structure, cutoffs, verbosity=4)

# orbitlist.sort()

for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, orbitlist.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(orbitlist.size()))

# for i in range(orbitlist.size()):
#      (orbitlist.get_orbit(i).get_representative_cluster().print())
# print("number of equivalent sites:
# ",len(orbitlist.get_orbit(i).get_equivalent_sites()))


#################################################
#                                               #
#        Get orbitlist for a supercell          #
#                                               #
#################################################

import random
N = 4
atoms = atoms.repeat(N)
for atom in atoms:
    if random.random() < 0.5:
        atom.symbol = "H"

structure_repeat = structure_from_atoms(atoms)
print("size of atoms {}".format(len(atoms)))
# for i in range(len(cutoffs) + 2):
#     print("number of {0}body clusters = {1}".format(
#         i, supercell_orbitlist.get_number_of_NClusters(i)))


#################################################
#                                               #
#             Count orbitlist                   #
#                                               #
#################################################

clustercounts1 = ClusterCounts()
t1 = time.time()
clustercounts1.count_each_local_orbitlist(
    structure_repeat, orbitlist_no_symmetry)
asd = clustercounts1.get_cluster_counts()
t2 = time.time()

print("Time to count each local orbitlist for orbitlist_no_symmetry {0:.6f} s".format(
    t2 - t1))

clustercounts1.print()

localClusterMap_no_symmetry = clustercounts1.get_cluster_counts()


clustercounts2 = ClusterCounts()
t1 = time.time()
clustercounts2.count_each_local_orbitlist(structure_repeat, orbitlist)
asd = clustercounts2.get_cluster_counts()
t2 = time.time()

print("Time to count each local orbitlist {0:.6f} s".format(t2 - t1))

clustercounts2.print()

localClusterMap = clustercounts2.get_cluster_counts()



t1 = time.time()
supercell_orbitlist = orbitlist.get_supercell_orbitlist(structure_repeat)
t2 = time.time()
supercell_orbitlist_time = t2 - t1
print("Time to get supercell orbitlist with {0} atoms : {1:.6f} s".format(
    structure_repeat.size(), t2 - t1))

print("size of repeated orbitlist {0}".format(supercell_orbitlist.size()))


clustercounts3 = ClusterCounts()
t1 = time.time()
clustercounts3.count_orbitlist(structure_repeat, supercell_orbitlist)
asd = clustercounts3.get_cluster_counts()
t2 = time.time()
print("Time to count orbitlist {0:.6f} s".format(t2 - t1))

print("Time to count orbitlist including setup {0:.6f} s".format(
    t2 - t1 + supercell_orbitlist_time))




allClusterMap = clustercounts3.get_cluster_counts()

for key in allClusterMap.keys():
    print(key.get_sites(), key.get_distances())
    print("")
    print("local")
    print(localClusterMap[key])
    print("no symmetry ol")
    print(localClusterMap_no_symmetry[key])

    print("global")
    print(allClusterMap[key])
    print("-----------")

for key in localClusterMap_no_symmetry.keys():
    print(key.get_sites(), key.get_distances())
    print("")

    print("no symmetry ol")
    print(localClusterMap_no_symmetry[key])
    print("-----------")

