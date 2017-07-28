from icetdev.orbitList import create_orbit_list
from ase import Atoms
from ase.build import bulk
from icetdev.permutationMap import PermutationMap, permutation_maps_from_atoms
from icetdev.structure import structure_from_atoms
from icetdev.manybodyNeighborlist import get_all_lattice_neighbors, ManybodyNeighborlist
import numpy as np
from icetdev.clusterCounts import ClusterCounts

from icetdev.orbitList import create_orbit_list, __get_latNbr_permutation_matrix
from time import time
atoms = bulk("Al", "fcc", a=2.0).repeat(1)


cutoffs = [10.5, 10]
pm_maps, prim_structure, neighborlists = permutation_maps_from_atoms(atoms, cutoffs, verbosity=0)



pm_matrix = pm_maps[0]

row1 = 1
row2 = 2
pm = pm_matrix.get_permutated_positions()
row1 = pm[row1]
row2 = pm[row2]

# for l1,l2 in zip(row1, row2):
#     print("dists : ",   np.linalg.norm(np.dot(l1-l2,prim_structure.cell)))
# # exit(1)

pm_lat_nbr = __get_latNbr_permutation_matrix(prim_structure, pm_matrix)


row1 = 1
row2 = 3

row1 = pm_lat_nbr[row1]
row2 = pm_lat_nbr[row2]

# for l1,l2 in zip(row1, row2):
#     # print(l1,l2)
#     print("dists : ",   prim_structure.get_distance2(l1.index, l1.unitcellOffset, l2.index, l2.unitcellOffset) )

# exit(1)    
# print("col1")
# for row in pm_lat_nbr:
#     print(row[0])
# print("col1 end")

# print("nl:")
# for latNbr in neighborlists[0].get_neighbors(0):
#     print(latNbr)
# print("nl end.")

clusters = create_orbit_list(structure=prim_structure, permutation_matrix=pm_matrix, neighborlists=neighborlists)
print("len of clusters ",len(clusters))
# #get lattice neighbors


clusterCounts = ClusterCounts()

clusterCounts.count_clusters(atoms=atoms, cutoffs=cutoffs)


print("Found {} clusters".format(clusterCounts.size()))



#lattice_neighbors = get_all_lattice_neighbors(structure=prim_structure, cutoffs=cutoffs)

# nbr_total = {}
# for i in range(5):
#     nbr_total[i] = 0
# print("len of latnbrs {}".format(len(lattice_neighbors)))
# for latNbrs in lattice_neighbors:    
#     try:
#         print(len(latNbrs[0]), len(latNbrs[1]))
#         nbr_total[len(latNbrs[0])] += len(latNbrs[1])
#     except:
#         print(len(latNbrs[1]))   


# print(nbr_total)             