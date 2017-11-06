from _icetdev import ClusterCounts
from icetdev.structure import structure_from_atoms
from icetdev.neighborlist import Neighborlist
from icetdev.manybody_neighborlist import ManybodyNeighborlist
from icetdev.cluster import Cluster
from icetdev.local_orbitlist_generator import LocalOrbitlistGenerator


def __count_each_local_orbitlist(self, structure, prim_orbitlist, keepOrderIntact=False):
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
    keepOrderIntact: bool
        count the clusters in the orbit with the same orientation as the prototype cluster    
    """

    localOrbitListGenerator = LocalOrbitlistGenerator(
        prim_orbitlist, structure)

    for i in range(localOrbitListGenerator.get_unique_offsets_count()):
        # sending local ol directly into function was about 10% faster
        # local_orbitlist = localOrbitListGenerator.generate_local_orbitlist(i)
        self.count_orbitlist(
            structure, localOrbitListGenerator.generate_local_orbitlist(i), keepOrderIntact)


ClusterCounts.count_each_local_orbitlist = __count_each_local_orbitlist
