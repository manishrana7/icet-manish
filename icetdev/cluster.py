from _icetdev import Cluster


def __len_of_cluster(self):
    return self.get_number_of_bodies()

Cluster.__len__ = __len_of_cluster
