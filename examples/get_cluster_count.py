from icetdev.clusterCounts import ClusterCounts
from ase.build import bulk
import numpy as np
atoms = bulk("Al").repeat(5)

atoms.set_chemical_symbols([['Al','Cu'][n] for n in np.round(np.random.random((len(atoms),))).astype(int)])

clusterCounts = ClusterCounts()

clusterCounts.count_clusters(atoms=atoms,cutoffs=2*[6])


clusterCounts.print()
#cluster_count_map = clusterCounts.get_cluster_counts()

#for cluster in cluster_count_map:
#    print(cluster)