"""
This scripts shows how to get print some basic information about the clusterspace
"""

from icetdev import clusterspace
from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms
from ase.build import bulk, make_supercell
import numpy as np

cutoffs = [10.0, 5.0]
subelements = ['Re', 'Ti']
prototype = bulk('Re')
clusterspace = create_clusterspace(subelements, cutoffs, atoms=prototype)

conf = structure_from_atoms(prototype.repeat(10))

cv = clusterspace.get_clustervector(conf)

print(cv)
print(clusterspace)