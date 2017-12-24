'''
This example demonstrates how to find the native clusters for a structure
'''
# Import modules
import numpy as np

from ase.build import bulk

from icetdev.cluster_space import ClusterSpace
from icetdev.structure import Structure


# Create a prototype structure, decide which additional elements to populate
# it with (Si, Ge) and set the cutoff for pairs (10.0 A)
conf = bulk("Si")
cutoffs = [10.0]
subelements = ["Si", "Ge"]

# Initiate the cluster space.
clusterspace = ClusterSpace(conf, cutoffs, subelements)

# Prepare 2x2x1 supercells, populate these, randomly, with Si and Ge atoms.
supercell = bulk("Si").repeat([2, 2, 1])
for atom in supercell:
    atom.symbol = np.random.choice(subelements)
structure = Structure.from_atoms(supercell)

# Extract and print the native clusters for the supercell.
nativeclusters = clusterspace.get_native_clusters(structure)
print(structure)
print("\nNative cluster counts:")
nativeclusters.print()
