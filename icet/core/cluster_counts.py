from ase import Atoms
from icet.core.orbit_list import OrbitList
from icet import Structure
from _icet import ClusterCounts as _ClusterCounts
from _icet import Cluster
from .local_orbit_list_generator import LocalOrbitListGenerator
from collections import OrderedDict


class ClusterCounts(_ClusterCounts):
    """
    Provides an interface to inspect cluster counts.

    Parameters
    ----------
    orbit_list : OrbitList object
        orbit list for a primitive structure
    atoms : ASE Atoms object
        supercell of the structure that `orbit_list` is based on

    Attributes
    ----------
    cluster_counts : OrderedDict
        keys are representative clusters (icet Cluster objects) for all
        distinct orbits in the orbit list, values are dicts where keys are
        the elements in a cluster, and values the number of counts of such
        clusters, e.g. {('Au', 'Ag'): 3, ('Au', 'Au'): 5}
    """

    def __init__(self, orbit_list: OrbitList, atoms: Atoms):
        self._orbit_list = orbit_list
        self._structure = Structure.from_atoms(atoms)

        # call (base) C++ constructor
        _ClusterCounts.__init__(self)

        self.cluster_counts = self._count_clusters()

    def _count_clusters(self, keep_order_intact=False, permute_sites=True):
        """
        Count all clusters in a structure by finding their local orbit list.

        Parameters
        ----------
        keep_order_intact: boolean
            if true the order in the cluster will be sorted
        permute_sites : boolean
            if true will permute the sites so they are in the
            symmetrically equivalent order as the representative sites
        """

        local_orbit_list_generator = LocalOrbitListGenerator(
            self._orbit_list, self._structure)

        for i in range(local_orbit_list_generator.get_unique_offsets_count()):
            self.count_orbit_list(
                self._structure,
                local_orbit_list_generator.generate_local_orbit_list(i),
                keep_order_intact, permute_sites)

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

        # Put clusters in a sorted order in an OrderedDict together with their
        # actual counts
        sorted_cluster_counts = OrderedDict()
        for cluster_data in sorted(cluster_counts):
            cluster = cluster_data[2]
            sorted_cluster_counts[cluster] = {}
            for m in range(*cluster_data[3]):
                elements, count = self.get_cluster_counts_info(m)
                sorted_cluster_counts[cluster][tuple(elements)] = count
        return sorted_cluster_counts

    def __str__(self):
        """
        String representation of cluster counts that provides an overview
        of the clusters (cluster, elements and count).
        """
        tuplets = {1: 'Singlet', 2: 'Pair', 3: 'Triplet', 4: 'Quadruplet'}
        width = 60
        s = ['{s:=^{n}}'.format(s=' Cluster Counts ', n=width)]

        first = True
        for cluster, counts in self.cluster_counts.items():
            if not first:
                s += ['']
            else:
                first = False

            # Add a description of the orbit to the string
            tuplet_type = tuplets.get(len(cluster.sites),
                                      '{}-tuplet'.format(len(cluster.sites)))
            s += ['{}: {} {:} {:.4f}'.format(tuplet_type,
                                             cluster.sites,
                                             cluster.distances,
                                             cluster.radius)]

            # Print the actual counts together with the elements they refer to
            for elements, count in counts.items():
                t = ['{:3} '.format(el) for el in elements]
                s += ['{} {}'.format(''.join(t), count)]
        s += [''.center(width, '=')]
        return '\n'.join(s)

    def __getitem__(self, key):
        """
        Returns cluster count (Cluster object and dict with counts) for a
        ClusterCounts object.

        Parameters
        ----------
        key : int / icet Cluster object (bi-optional)
            If int, return the key:th counts, if Cluster, return the counts
            belonging to that cluster.
        """
        if isinstance(key, int):
            return list(self.cluster_counts.values())[key]
        elif isinstance(key, Cluster):
            return self.cluster_counts[key]
