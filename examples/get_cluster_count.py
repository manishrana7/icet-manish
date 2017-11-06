from icetdev.cluster_counts import ClusterCounts
from ase.build import bulk
import numpy as np
atoms = bulk("Ti","bcc", a=3.43).repeat([2,2,1])
atoms.pbc = [True,True,False]
atoms.set_chemical_symbols([['Ti','W'][n] for n in np.round(np.random.random((len(atoms),))).astype(int)])

cluster_counts = ClusterCounts()
cutoffs = [4]
cluster_counts.count_clusters(atoms=atoms, cutoffs=cutoffs)

print("number of atoms {0}".format(len(atoms)))
print("Found {} clusters".format(cluster_counts.size()))
cluster_counts.print()
#cluster_count_map = cluster_counts.get_cluster_counts()

#for cluster in cluster_count_map:
#    print(cluster)
