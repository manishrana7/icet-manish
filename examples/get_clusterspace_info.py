"""
This scripts shows how to get print some basic information about the clusterspace
"""

from icetdev import clusterspace
from icetdev.clusterspace import create_clusterspace, get_singlet_info, view_singlets
from icetdev.structure import structure_from_atoms
from ase.build import bulk, make_supercell
import numpy as np

cutoffs = [10.0, 7.0]

subelements = ['Re', 'Ti',"H"]
prototype = bulk('Re',"bcc",a=3).repeat(3)
prototype.pbc = [False, False, False]
clusterspace = create_clusterspace(prototype, cutoffs, subelements)

#Get clusterspace details by printing
print(clusterspace)