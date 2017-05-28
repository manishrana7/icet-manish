import icetdev as it
import numpy as np

sites = np.array([0, 0, 0])
distances = np.array([3.0, 1.0, 0.0])

cluster = it.Cluster(sites, distances)


assert cluster.get_count([0, 0, 0]) == 0, "if no count, then its zero"

#count
cluster.count([0, 0, 0])

assert cluster.get_count([0, 0, 0]) == 1, "first count is one"


assert (cluster.get_distances() == distances).all(), "The distances should be equal"

assert (cluster.get_sites() == sites).all(), "The sites should be equal"


#test equality

sites_big = np.array([10,10,10])
distances_big = np.array([10,10,10])
cluster_big = it.Cluster(sites_big, distances_big)

assert cluster < cluster_big

#try distances equal but sites different

sites_slightly_bigger = np.array([11,10,10])
distances_slightly_bigger = np.array([10,10,10])
cluster_slightly_bigger = it.Cluster(sites_slightly_bigger, distances_slightly_bigger)


assert cluster_big < cluster_slightly_bigger

#try distances different but sites same

sites_slightly_bigger = np.array([10,10,10])
distances_slightly_bigger = np.array([10.1,10,10])
cluster_slightly_bigger = it.Cluster(sites_slightly_bigger, distances_slightly_bigger)

assert cluster_big < cluster_slightly_bigger

#make sure that the same cluster is not less then
cluster_equal = it.Cluster(sites_slightly_bigger, distances_slightly_bigger)
assert not (cluster_slightly_bigger < cluster_equal)

assert not (cluster_equal < cluster_slightly_bigger)

#make sure that smaller bodied cluster are consider smaller than higher bodied clusters

sites_few = np.array([100,1000])
distances_few = np.array([1000000.0])
cluster_pair = it.Cluster(sites_few, distances_few)
assert cluster_pair < cluster_big