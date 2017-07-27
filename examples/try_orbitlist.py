from icetdev.orbitList import create_orbit_list
from ase import Atoms
from ase.build import bulk
from icetdev.permutationMap import PermutationMap, permutation_maps_from_atoms
from icetdev.structure import structure_from_atoms
import numpy as np


atoms = bulk("Al", "fcc", a=2.0).repeat(1)


cutoffs = [10]
pm_maps, prim_structure = permutation_maps_from_atoms(atoms, cutoffs, verbosity=0)

pm_matrix = pm_maps[0]

#get mbnl....

