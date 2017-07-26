from icetdev.tools.geometry import get_scaled_positions, find_latticeNeighbor_from_position

from icetdev.structure import *
from icetdev.neighborlist import *
from ase import Atoms, Atom
from ase.build import bulk

import numpy as np


# ASE atoms
atoms = bulk("Al", "fcc", a=2.0).repeat(1)
atoms.append(Atom("Al",position=[0.5,0.5,0.5]))

#get neighborlist
neighbor_cutoff = 4.0 #Ångström
structure = structure_from_atoms(atoms)
nl = Neighborlist(neighbor_cutoff)
nl.build(structure)
pos_neighbors = []

#get frac positions
for latNbr in nl.get_neighbors(0):
    pos = structure.get_position(latNbr)
    pos_neighbors.append(pos)
frac_coordinates = get_scaled_positions(np.array(pos_neighbors), cell=atoms.cell, wrap=False, pbc = structure.pbc)


#print to see what happens
print("#1: frac pos #2: position #3: fpos*cell")
for fpos, pos in zip(frac_coordinates,pos_neighbors):
    lat_nbr = find_latticeNeighbor_from_position(structure, pos)
    lat_pos = structure.positions[lat_nbr.index] + np.dot(lat_nbr.unitcellOffset, structure.cell)
    print(fpos, pos, np.dot(fpos, atoms.cell), lat_nbr, lat_pos-pos )


