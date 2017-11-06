from icetdev.clusterspace import create_clusterspace
from icetdev.cluster_counts import ClusterCounts
from icetdev.structure import structure_from_atoms
from ase.build import bulk
from ase import Atoms
import numpy as np




def print_native_clusters(clustercounts, structure):
    print("Native cluster counts for:")
    print(structure)
    clustercounts.print()
    

conf = bulk("Si")
cutoffs = [10.0]
subelements = ["Si", "Ge"]

clusterspace = create_clusterspace(conf, cutoffs, subelements)


for i in range(1):
    supercell = bulk("Si").repeat([2,2,1])
    for atom in supercell:
        atom.symbol = np.random.choice(subelements)
    supercell = structure_from_atoms(supercell)        
    clustercounts = clusterspace.get_native_clusters(supercell)
    print_native_clusters(clustercounts, supercell)
