from ase import Atoms
from icet.core.orbit_list import OrbitList
from icet import Structure
from _icet import ClusterCounts as _ClusterCounts
from .local_orbit_list_generator import LocalOrbitListGenerator


class ClusterCounts(_ClusterCounts):
    """
    Provides an interface to inspect cluster counts.

    Parameters
    ----------
    primitive_orbit_list : OrbitList object
        orbit list for a primitive structure
    atoms : ASE Atoms object
        supercell of the primitive structure that `primitive_orbit_list` is
        based on
    """

    def __init__(self, primitive_orbit_list: OrbitList, atoms: Atoms):
        self._primitive_orbit_list = primitive_orbit_list
        self._structure = Structure.from_atoms(atoms)

        # call (base) C++ constructor
        _ClusterCounts.__init__(self)

        self.cluster_counts = self._count_clusters()

    def _count_clusters(self, keep_order_intact=False):
        """
        Count all clusters in a structure by finding their local orbit list.

        Parameters
        ----------
        keep_order_intact: boolean
            count the clusters in the orbit with the same orientation as the
            prototype cluster
        """
        local_orbit_list_generator = LocalOrbitListGenerator(
            self._primitive_orbit_list, self._structure)

        for i in range(local_orbit_list_generator.get_unique_offsets_count()):
            # sending local orbit list directly into function was about
            # 10% faster than:
            # local_orbit_list = \
            #    local_orbit_list_generator.generate_local_orbit_list(i)
            self.count_orbit_list(
                self._structure,
                local_orbit_list_generator.generate_local_orbit_list(i),
                keep_order_intact)

        self.setup_cluster_counts_info()

        # Put clusters in a list that we can sort according to
        # (1) multiplet and (2) size of cluster
        m = 0
        cluster_counts = []
        for cluster, cluster_info in self.get_cluster_counts().items():
            cluster_counts.append(
                (len(cluster.distances), cluster.radius, cluster,
                 (m, m+len(cluster_info))))
            m += len(cluster_info)

        # Put clusters in a sorted order in a list together with their
        # actual counts
        sorted_cluster_counts = []
        for cluster_data in sorted(cluster_counts):
            sorted_cluster_counts.append((cluster_data[2], {}))
            for m in range(*cluster_data[3]):
                elements, count = self.get_cluster_counts_info(m)
                sorted_cluster_counts[-1][1][tuple(elements)] = count
        return sorted_cluster_counts

    def __repr__(self):
        """
        String representation of cluster counts that provides an overview
        of the clusters (cluster, elements and count).
        """
        tuplets = {1: 'Singlet', 2: 'Pair', 3: 'Triplet', 4: 'Quadruplet'}
        width = 60
        s = ['{s:=^{n}}'.format(s=' Cluster Counts ', n=width)]

        first = True
        for cluster_count in self.cluster_counts:
            if not first:
                s += ['']
            else:
                first = False

            cluster = cluster_count[0]

            # Add a description of the orbit to the string
            tuplet_type = tuplets.get(len(cluster.sites),
                                      '{}-tuplet'.format(len(cluster.sites)))
            s += ['{}: {} {:} {:.4f}'.format(tuplet_type,
                                             cluster.sites,
                                             cluster.distances,
                                             cluster.radius)]

            # Print the actual counts together with the elements they refer to
            for elements, count in cluster_count[1].items():
                t = ['{:3} '.format(el) for el in elements]
                s += ['{} {}'.format(''.join(t), count)]
        s += [''.center(width, '=')]
        return '\n'.join(s)

    def __getitem__(self, key):
        """
        Return cluster count (Cluster object and dict with counts) for a
        ClusterCounts object.
        """
        return self.cluster_counts[key]
