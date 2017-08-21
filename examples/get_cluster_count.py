from icetdev.clusterCounts import ClusterCounts
from ase.build import bulk
import numpy as np
atoms = bulk("Ti","bcc", a=3.43).repeat([2,2,1])
atoms.pbc = [True,True,False]
atoms.set_chemical_symbols([['Ti','W'][n] for n in np.round(np.random.random((len(atoms),))).astype(int)])

clusterCounts = ClusterCounts()
cutoffs = [4]
clusterCounts.count_clusters(atoms=atoms, cutoffs=cutoffs)

print("number of atoms {0}".format(len(atoms)))
print("Found {} clusters".format(clusterCounts.size()))
clusterCounts.print()
#cluster_count_map = clusterCounts.get_cluster_counts()

#for cluster in cluster_count_map:
#    print(cluster)