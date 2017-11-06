from _icetdev import ClusterCounts
from icetdev.local_orbitlist_generator import LocalOrbitlistGenerator


def __count_clusters(self, structure,
                     prim_orbitlist, keepOrderIntact=False):
    """
    Count all clusters in a structure by finding their local orbitlist.

    This will need to create the localOrbitListGenerator object for generating
    the local orbitlists.

    Parameters
    ----------
    structure : icet structure object
        supercell of the structure prim_orbitlist is based on
    prim_orbitlist : icet orbitlist object
        based on a primitive of the input structure
    keepOrderIntact: bool
        count the clusters in the orbit with the same orientation as the
        prototype cluster
    """

    localOrbitListGenerator = LocalOrbitlistGenerator(
        prim_orbitlist, structure)

    for i in range(localOrbitListGenerator.get_unique_offsets_count()):
        # sending local ol directly into function was about 10% faster
        # local_orbitlist = localOrbitListGenerator.generate_local_orbitlist(i)
        self.count_orbitlist(
            structure,
            localOrbitListGenerator.generate_local_orbitlist(i),
            keepOrderIntact)


ClusterCounts.count_clusters = __count_clusters
