'''
This example demonstrates how to find the native clusters for a structure
'''

# Import modules
import numpy as np
from ase.build import bulk
from icet import ClusterSpace, Structure

# Create a prototype structure, decide which additional elements to populate
# it with (Si, Ge) and set the cutoff for pairs (10.0 A)
conf = bulk('Si')
cutoffs = [10.0]
subelements = ['Si', 'Ge']

# Initiate the cluster space.
cluster_space = ClusterSpace(conf, cutoffs, subelements)

# Prepare 2x2x1 supercells, populate these, randomly, with Si and Ge atoms.
supercell = bulk('Si').repeat([2, 2, 1])
for atom in supercell:
    atom.symbol = np.random.choice(subelements)

# Extract and print the native clusters for the supercell.
# TODO: The cluster_space.get_native_clusters method is provided
# directly via Pybind. What does it do? Is it necessary to expose it
# to the user? This appears to be the _only_ place in Python where it
# shows up.
structure = Structure.from_atoms(supercell)
native_clusters = cluster_space.get_native_clusters(structure)
print(structure)
print('\nNative cluster counts:')
native_clusters.print()
