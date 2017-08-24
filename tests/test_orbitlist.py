from icetdev import orbitList
from icetdev.structure import structure_from_atoms, Structure
from icetdev.manybodyNeighborlist import *
from ase import Atoms
from ase.build import bulk
from icetdev.neighborlist import get_neighborlists
from icetdev.orbitList import create_orbit_list


def setup_test_orbitlist(atoms, cutoffs):
    structure = structure_from_atoms(atoms)
    mbnl = ManybodyNeighborlist()
    neighborlists = get_neighborlists(atoms=atoms, cutoffs=cutoffs)
    mbnl.build(neighborlists, 0, True)
    return structure, mbnl


atoms = bulk("Al", "diamond", a=1)
cutoffs = [3.01, 2.01]


structure, mbnl = setup_test_orbitlist(atoms, cutoffs)

ol = orbitList.OrbitList(mbnl, structure)

ol.sort()

#ol.print(verbosity = 3)
for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, ol.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(ol.size()))


############################################
#                                          #
#   Test orbitlist from permutation map    #
#                                          #
############################################

orbitlist = create_orbit_list(structure, cutoffs, verbosity=4)

# orbitlist.sort()

#ol.print(verbosity = 3)
for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, orbitlist.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(orbitlist.size()))

import numpy as np
for i in range(orbitlist.size()):
    if len(orbitlist.get_orbit(i).get_representative_cluster().get_distances())==3:        
        print(orbitlist.get_orbit(i).get_representative_cluster().get_distances())


#         eq_sites = orbitlist.get_orbit(i).get_equivalent_sites()
#         sites = eq_sites[0]
#         for s in sites:
#             print(np.dot(s.unitcellOffset, atoms.cell))
#         print(" end ")
#         print(" ")
# for i in range(orbitlist.size()):
#     if len(orbitlist.get_orbit(i).get_representative_cluster().get_distances())==3:        
#         dists = orbitlist.get_orbit(i).get_representative_cluster().get_distances()
#         if np.linalg.norm((np.array(dists) -[0.86, 1.4,1.65])) < 0.1:
#             eq_sites = orbitlist.get_orbit(i).get_equivalent_sites()
#             print("len of eq sites: {}".format(len(eq_sites)) )
#             for sites in eq_sites:
#                 for s in sites:                    
#                     print(np.dot(s.unitcellOffset, atoms.cell), end= ' ')
#                 print( " ")                    
#             print(" end ")
#             print(" ")
