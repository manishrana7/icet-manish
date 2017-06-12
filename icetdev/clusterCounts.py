from example import ClusterCounts
from icetdev.structure import structure_from_atoms
from icetdev.neighborlist import Neighborlist
from icetdev.manybodyNeighborlist import ManybodyNeighborlist


def __count_clusters(self, atoms=None, structure=None, nl=None, mbnl=None, cutoffs=None,
                     reset=True, order=None):
    """
    Counts clusters in either atoms or structure

    Parameters
    ----------
    atoms:
        ASE atoms object (bi-optional) need atleast this or structure
    structure:
        icet structure object (bi-optional) need atleast this or atoms
    nl:
        icet neighborlist object
        if None: then it will be created
    mbnl:
        icet manybody neighborlist
        if None: then it will be created
    cutoffs:
        positive floats indicating the cutoffs for the various clusters
        it is not used if nl != None
        this will set the order if order is None
        order will be set to: len(cutoffs) + 1
        (plus one since the first cutoff is for pairs)
    reset:
        boolean : default = True
        will reset clustercount if true
    order:
        int : needed if cutoffs is not set
        sets the maximum cluster order
    """

    bothways = False
    if reset:
        self.reset()
    if structure is None:
        if atoms is None:
            raise Exception(
                "Error: both atoms and structure is None in count clusters")
        else:
            structure = structure_from_atoms(atoms)

    if nl is None:
        if cutoffs is None:
            raise Exception(
                "Error: both nl and cutoffs is None in count clusters")
        else:
            nl = Neighborlist(max(cutoffs))
    nl.build(structure)
    if cutoffs is None:
        if order is None:
            raise Exception(
                "Error: both order and cutoffs is None in count clusters")
    else:
        order = len(cutoffs) + 1
    if mbnl is None:
        mbnl = ManybodyNeighborlist()

    # loop over the indices in the structure, create the mbnl for each index
    # and count each index
    for lattice_index in range(structure.size()):
        lattice_neigbhors = mbnl.build(nl, lattice_index, order, bothways)
        # for latNbrs in lattice_neigbhors:
        #     for latnbrComb in latNbrs[1]:
        #         for latnbr in latNbrs[0]:
        #             print(latnbr.index, latnbr.unitcellOffset, end=' ')
        #         print(latnbrComb.index, latnbrComb.unitcellOffset)
        self.count_lattice_neighbors(structure, lattice_neigbhors)
    # return all objects that might have been created here for possible reuse
    return structure, nl, mbnl, order


ClusterCounts.count_clusters = __count_clusters
