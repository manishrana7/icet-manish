import numpy as np
from ase.db import connect
from ase.build import bulk, cut
from icetdev.cluster_space import ClusterSpace
from ase.visualize import view
from ase.io import write
from spglib import get_spacegroup
from spglib import niggli_reduce
from icetdev.permutation_map import __get_primitive_structure
from icetdev.tools.geometry import transform_cell_to_cell, required_offsets_to_map_supercell, get_smart_offsets

from icetdev.tools.map_sites import map_configuration_to_reference



p_trial=[[1, 0, 0], [0, 1, 100], [0, 0, 2]]

# initialize cluster space and get the internal primitive atoms object
prim=bulk('Au', a=4.0, crystalstructure='hcp')
subelements=['Au', 'Pd']
cutoffs=[8.0, 6.0]
cs=ClusterSpace(prim, cutoffs, subelements)

atoms_prim=cs.get_primitive_structure().to_atoms()


# Create a supercell using permutation matrix
atoms=cut(atoms_prim, p_trial[0], p_trial[1], p_trial[2])
print("# Len of primitive: ",len(atoms_prim))
print("# Len of supercell: ",len(atoms))


unique_offsets = required_offsets_to_map_supercell(atoms, atoms_prim)


print("The unique primitive cell offsets the supercell atoms end up in:")
print(unique_offsets)
print()
print("Will offsetting the primitive atoms object with these offsets create too many atoms? \nAnswer: {}\n ".format(len(unique_offsets) > len(atoms)/len(atoms_prim)))

print("Trying to find {} number of mappings of the primitive that maps the supercell once and only once...".format(len(atoms)//len(atoms_prim)))

smart_offsets = get_smart_offsets(atoms, atoms_prim)

print(smart_offsets)