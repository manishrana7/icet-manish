from icetdev.tools.geometry import get_scaled_positions, find_latticeNeighbor_from_position_python

from icetdev.structure import *
from icetdev.neighborlist import *
from ase import Atoms, Atom
from ase.build import bulk

import numpy as np


# ASE atoms
atoms = bulk("Al", "fcc", a=2.0).repeat(1)
atoms.append(Atom("Al",position=[0.5,0.5,0.5]))
atoms.append(Atom("Al",position=[0.25,0.5,0.5]))

#get neighborlist
neighbor_cutoff = 6.0 #Ångström
structure = structure_from_atoms(atoms)
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
pos_neighbors = []

#get frac positions
for latNbr in nl.get_neighbors(0):
    pos = structure.get_position(latNbr)
    pos_neighbors.append(pos)
frac_coordinates = get_scaled_positions(np.array(pos_neighbors), cell=atoms.cell, wrap=False, pbc = structure.pbc)

import time

t0 = time.time()
#print to see what happens
print("#1: frac pos #2: position #3: fpos*cell")
lat_nbrs = structure.findLatticeNeighborsFromPositions(pos_neighbors)

# for fpos, pos in zip(frac_coordinates,pos_neighbors):
#     lat_nbr = structure.findLatticeNeighborFromPosition(pos)
#     #lat_nbr = find_latticeNeighbor_from_position_python(structure, pos)
#     #lat_pos = structure.positions[lat_nbr.index] + np.dot(lat_nbr.unitcellOffset, structure.cell)
#     #print(fpos, pos, np.dot(fpos, atoms.cell), lat_nbr, lat_pos-pos )

# for fpos, pos, lat_nbr in zip(frac_coordinates,pos_neighbors, lat_nbrs):
#     #lat_nbr = structure.findLatticeNeighborFromPosition(pos)
#     #lat_nbr = find_latticeNeighbor_from_position_python(structure, pos)
#     lat_pos = structure.positions[lat_nbr.index] + np.dot(lat_nbr.unitcellOffset, structure.cell)
#     print(fpos, pos, np.dot(fpos, atoms.cell), lat_nbr, lat_pos-pos )



t1 = time.time()
print("time for finding latnbrs: ", round((t0-t1)*1e3,4))
print("Found {} LatticeNeighbors".format(len(lat_nbrs)))