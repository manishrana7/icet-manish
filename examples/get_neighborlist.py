import numpy as np
import numpy.random as random

from ase import Atoms
from ase.build import bulk
from ase.neighborlist import NeighborList
from icetdev import *
from icetdev.structure import *

neighbor_cutoff = 1.5

# ASE neighborlist
atoms = bulk("Al", "fcc", a=2).repeat(2)
atoms.pbc = [True, True, True]
structure = structure_from_atoms(atoms)


# icet neighborlist
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)

# neighbors come in a list of tuples
for index in range(len(atoms)):
    # get neighbors of index "index"
    neighbors = nl.get_neighbors(index)
    print("Neighbors of atom with index {}".format(index))
    for neighbor in neighbors:
        neighbor_index = neighbor[0]
        neighbor_offset = neighbor[1]
        distance_to_neighbor = structure.get_distance2(
            index, [0, 0, 0], neighbor[0], neighbor[1])
        print("{0} {1} {2:1.5f}".format(neighbor_index,
                                        neighbor_offset, distance_to_neighbor))
    print("")
print("fcc has {} nearest neighbors".format(len(neighbors)))
