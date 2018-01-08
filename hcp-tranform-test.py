import numpy as np
from ase.db import connect
from ase.build import bulk, cut
from icetdev import ClusterSpace
from ase.visualize import view
from ase.io import write
from spglib import get_spacegroup
from spglib import niggli_reduce
from icetdev.tools.geometry import transform_cell_to_cell, required_offsets_to_map_supercell, get_smart_offsets


p_trial = [[1, 0, 0], [3, 1, 20], [0, 0, 2]]

# initialize cluster space and get the internal primitive atoms object
prim = bulk('Au', a=4.0, crystalstructure='hcp')
subelements = ['Au', 'Pd']
cutoffs = [8.0, 6.0]
cs = ClusterSpace(prim, cutoffs, subelements)

atoms_prim = cs.get_primitive_structure().to_atoms()


# Create a supercell using permutation matrix
atoms = cut(atoms_prim, p_trial[0], p_trial[1], p_trial[2])
print("# Len of primitive: ", len(atoms_prim))
print("# Len of supercell: ", len(atoms))


unique_offsets = required_offsets_to_map_supercell(atoms, atoms_prim)


print("The unique primitive cell offsets the supercell atoms end up in:")
print(unique_offsets)
print()
print("Will offsetting the primitive atoms object with these offsets create too many atoms? \nAnswer: {}\n ".format(
    len(unique_offsets) > len(atoms) / len(atoms_prim)))

print("Trying to find {} number of mappings of the primitive that maps the supercell once and only once...".format(
    len(atoms) // len(atoms_prim)))

smart_offsets = get_smart_offsets(atoms, atoms_prim)
print("Smart offsets found:")
for so in smart_offsets:
    for i, j, offset in zip(so.prim_indices, so.super_indices, so.offsets):
        print("Primitive index ", i, " via unitcell offset ",
              offset, "resulted in supercell index: ", j)
        prim_translated_pos = atoms_prim[i].position + \
            np.dot(offset, atoms_prim.cell)
        diff = prim_translated_pos - atoms[j].position
        assert np.linalg.norm(
            diff) < 1e-4, " Atom position was not translated correctly: {}".format(diff)

print("Did we get correct number of offsets? ", len(
    atoms) // len(atoms_prim) == len(smart_offsets))

print("Did all supercell atoms get mapped?", end=' ')

supercell_indices_mapped = []
for so in smart_offsets:
    for indices in so.super_indices:
        supercell_indices_mapped.append(indices)
print(len(supercell_indices_mapped) == len(atoms) and len(
    set(supercell_indices_mapped)) == len(supercell_indices_mapped))
