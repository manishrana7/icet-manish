from _icetdev import Cluster



def __len_of_cluster(self):
    """Length of cluster a.k.a. number of bodies"""

    return self.get_number_of_bodies()

Cluster.__len__ = __len_of_cluster
