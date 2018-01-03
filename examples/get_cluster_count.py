'''
This example demonstrates how to count the number of clusters for a structure.
'''

# Start import
from ase.build import bulk
from icetdev import Structure
from icetdev.tools import get_primitive_structure
from icetdev.core.cluster_counts import ClusterCounts
from icetdev.core.orbit_list import create_orbit_list
# End import

# Create a titanium, single-layered, sheet and randomly populate some of the
# sites with W atoms.
# Start setup
atoms = bulk('Ti', 'bcc', a=3.43).repeat([2, 2, 1])
atoms.pbc = [True, True, False]
atoms.set_chemical_symbols(['Ti', 'W', 'W', 'Ti'])
# End setup

# Determine the orbitlist for the corresponding primitive structure for all
# pair clusters withhin the cutoff distance (4 angstrom).
# Start orbitlist
cutoffs = [4]
prim_atoms = get_primitive_structure(atoms)
prim_structure = Structure.from_atoms(prim_atoms)
prim_orbitlist = create_orbit_list(prim_structure, cutoffs)
# End orbitlist

# Use the primitive orbitlist to count the number of clusters.
# Start counting
structure = Structure.from_atoms(atoms)
cluster_counts = ClusterCounts()
cluster_counts.count_clusters(structure, prim_orbitlist)
# End counting

# Print all of the clusters that were found.
# Start results
print("number of atoms {0}".format(len(atoms)))
print("Found {} clusters".format(len(cluster_counts)))
cluster_counts.print()
# End results
