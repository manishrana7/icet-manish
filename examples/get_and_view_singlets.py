"""
This scripts shows how to get print some basic information about the clusterspace
"""

from icetdev import clusterspace
from icetdev.clusterspace import create_clusterspace, get_singlet_info, view_singlets
from icetdev.structure import structure_from_atoms
from ase.build import bulk, make_supercell
import numpy as np


#create a nanoparticle, i.e. pbc=false
prototype = bulk('Re',"bcc",a=3).repeat(3)
prototype.pbc = False

cluster_data = get_singlet_info(prototype)

#print singlet information
for singlet in cluster_data:
    for key in singlet.keys():
        print(key, singlet[key])
    print("")        


#Visually inspect singlets in this structure
view_singlets(prototype)