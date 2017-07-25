from icetdev.clusterCounts import ClusterCounts
from ase.build import bulk
import numpy as np
atoms = bulk("Ti",a=3.43).repeat(1)

atoms.set_chemical_symbols([['Ti','W'][n] for n in np.round(np.random.random((len(atoms),))).astype(int)])

clusterCounts = ClusterCounts()

clusterCounts.count_clusters(atoms=atoms, cutoffs=[15.0, 13,15,10])


print("Found {} clusters".format(clusterCounts.size()))
#clusterCounts.print()
#cluster_count_map = clusterCounts.get_cluster_counts()

#for cluster in cluster_count_map:
#    print(cluster)