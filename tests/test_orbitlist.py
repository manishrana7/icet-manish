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


atoms = bulk("Al", "bcc", a=1)
cutoffs = [2]


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

orbitlist = create_orbit_list(structure, cutoffs)

#orbitlist.sort()

#ol.print(verbosity = 3)
for i in range(len(cutoffs) + 2):
    print("number of {0}body clusters = {1}".format(
        i, orbitlist.get_number_of_NClusters(i)))

print("size of orbitlist {0}".format(orbitlist.size()))
