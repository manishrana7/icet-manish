class Clusters(object):
    def __init__(self, prototype, elements, cutoffs, 
                 basis_functions='standard'):
        self.prototype = prototype
        self.elements = elements
        self.cutoffs = cutoffs

        # obtain a cluster list using ClusterLib



    def get_number_of_distinct_clusters(self, ntuple=None):
        if ntuple == None:
            return len(self.cluster_list)
        else:
            return 0 # count only pairs, triplets etc

    def get_cluster_filter(self, cutoffs):
        return [0, 1, 2, 3, 4]

    def copy(self):
        return copy.copy(self)

    def copy_with_filter(self, cluster_filter):
        new_clusters = copy.copy(self)
        # make necessary rearrangements
        return new_clusters

    def save(self):
    	return 0

    def __len__(self):
        return len(self.cluster_list)

    def __str__(self):
        return('Clusters object with elements {} and '
               'cutoffs {}'.format(self.elements, self.cutoffs))
