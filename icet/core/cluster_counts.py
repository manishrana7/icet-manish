from _icet import ClusterCounts
from .local_orbit_list_generator import LocalOrbitListGenerator


def __count_clusters(self, structure, prim_orbit_list,
                     keep_order_intact=False):
    '''
    Count all clusters in a structure by finding their local orbit list.

    Note that this requires creating the local_orbit_list_generator object for
    generating the local orbit lists.

    Parameters
    ----------
    structure : Structure object
        supercell of the structure `prim_orbit_list` is based on
    prim_orbit_list : OrbitList object
        orbit list based on a primitive of the input structure
    keep_order_intact: boolean
        count the clusters in the orbit with the same orientation as the
        prototype cluster
    '''

    local_orbit_list_generator = LocalOrbitListGenerator(prim_orbit_list,
                                                         structure)

    for i in range(local_orbit_list_generator.get_unique_offsets_count()):
        # sending local orbit list directly into function was about 10% faster
        # than:
        # local_orbit_list = \
        #    local_orbit_list_generator.generate_local_orbit_list(i)
        self.count_orbit_list(
            structure,
            local_orbit_list_generator.generate_local_orbit_list(i),
            keep_order_intact)


def __get_string_representation(self):
    '''
    String representation of cluster counts that provides an overview
    of the clusters (cluster, elements and count).
    '''
    self.setup_cluster_counts_info()
    tuplets = {1: 'Singlet', 2: 'Pair', 3: 'Triplet', 4: 'Quadruplet'}
    width = 50
    s = ['{s:=^{n}}'.format(s=' Cluster Counts ', n=width)]
    cluster_counts = {key: len(values)
                      for key, values in self.get_cluster_counts().items()}
    m = 0
    horizontal_line = '{s:-^{n}}'.format(s='', n=width)
    for cluster in sorted(cluster_counts.keys()):
        if m != 0:
            s += ['']
        tuplet_type = tuplets.get(len(cluster.sites),
                                  '{}-tuplet'.format(len(cluster.sites)))
        s += ['{}: {} {:} {:.4f}'.format(tuplet_type,
                                         cluster.sites,
                                         cluster.distances,
                                         cluster.radius)]
        for k in range(cluster_counts[cluster]):
            elements, count = self.get_cluster_counts_info(m)
            t = ['{:3} '.format(el) for el in elements]
            s += ['{} {}'.format(''.join(t), count)]
            m += 1
    s += [''.center(width, '=')]
    return '\n'.join(s)


def __print_overview(self):
    '''
    Print cluster counts representation
    '''
    print(__get_string_representation(self))


ClusterCounts.count_clusters = __count_clusters
ClusterCounts.__repr__ = __get_string_representation
ClusterCounts.print_overview = __print_overview
