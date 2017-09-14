from _icetdev import ClusterCounts
from icetdev.structure import structure_from_atoms
from icetdev.neighborlist import Neighborlist
from icetdev.manybodyNeighborlist import ManybodyNeighborlist
from icetdev.cluster import Cluster
from icetdev.localOrbitlistGenerator import LocalOrbitlistGenerator


def __count_orbitlist(self, structure, orbitlist):
    """
    Counts the clusters of the sites in each orbit of orbitlist

    Parameters
    ----------
    structure:
        icet structure object to be counted on

    orbitlist:
        icet orbitlist class
    """

    for i, orbit in enumerate(orbitlist.get_orbitList()):
        repr_cluster = orbit.get_representative_cluster()
        cluster = Cluster(distances=repr_cluster.get_distances(),
                          sites=repr_cluster.get_sites(),
                          sortedCluster=False,
                          clusterTag=i)
        self.count(structure, orbit.get_equivalent_sites(), cluster)


ClusterCounts.count_orbitlist = __count_orbitlist


def __count_each_local_orbitlist(self, structure, prim_orbitlist):
    """
    Span the structure and count all clusters by finding each local orbitlist.

    This will need to create the localOrbitListGenerator object which can generate
    the local orbitlists

    Parameters
    ----------
    structure:
        icet structure object (supercell of the structure prim_orbitlist was based on)
    prim_orbitlist:        
        icet orbitlist object (based on a primitive of the input structure)
    """

    localOrbitListGenerator = LocalOrbitlistGenerator(
        prim_orbitlist, structure)

    for i in range(localOrbitListGenerator.get_unique_offsets_count()):
        # sending local ol directly into function was about 10% faster
        # local_orbitlist = localOrbitListGenerator.generate_local_orbitlist(i)
        self.count_orbitlist(
            structure, localOrbitListGenerator.generate_local_orbitlist(i))


ClusterCounts.count_each_local_orbitlist = __count_each_local_orbitlist
