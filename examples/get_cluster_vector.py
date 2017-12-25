'''
This example demonstrates how to construct cluster vectors.
'''
# Import modules
from ase.build import bulk

from icetdev.cluster_space import ClusterSpace

# Create a prototype structure, decide which additional elements to populate
# it with (Si, Ge) and set the cutoffs for pairs (5.0 Å), triplets (5.0 Å)
# and quadruplets (5.0 Å).
conf = bulk("Si")
cutoffs = [5.0, 5.0, 5.0]
subelements = ["Si", "Ge"]

# Initiate and print the cluster space.
clusterspace = ClusterSpace(conf, cutoffs, subelements)
print(clusterspace)

# Generate and print the cluster vector for a pure Si 2x2x2 supercell.
supercell = bulk("Si").repeat(2)
cv = clusterspace.get_cluster_vector(supercell)
print(cv)

# Generate and print the cluster vector for a mixed Si-Ge 2x2x2 supercell
supercell_2 = bulk("Si").repeat(2)
supercell_2[0].symbol = "Ge"
cv_2 = clusterspace.get_cluster_vector(supercell_2)
print(cv_2)
