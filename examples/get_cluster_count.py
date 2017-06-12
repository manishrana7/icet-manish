from icetdev.clusterCounts import ClusterCounts
from ase.build import bulk

atoms = bulk("Al").repeat(2)


clusterCounts = ClusterCounts()

clusterCounts.count_clusters(atoms=atoms,cutoffs=[5.0,5.0])


clusterCounts.print()
#cluster_count_map = clusterCounts.get_cluster_counts()

#for cluster in cluster_count_map:
#    print(cluster)