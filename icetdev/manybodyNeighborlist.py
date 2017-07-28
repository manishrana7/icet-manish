from example import ManybodyNeighborlist

from icetdev.latticeNeighbor import LatticeNeighbor
from icetdev.neighborlist import get_neighborlists


def get_all_lattice_neighbors(atoms=None, structure=None, neighborlists=None, cutoffs=None):
    """
     Return lattice neighbors from a configuration and cutoffs

    return list of lattice neighbors in mbnl format, list(neighborlists)
    Parameters
    ----------
    atoms:
        ASE atoms object (bi-optional) need atleast this or structure
    structure:
        icet structure object (bi-optional) need atleast this or atoms
    neighborlists:
        array/list of icet neighborlist object
        if None: then it will be created
    cutoffs:
        positive floats indicating the cutoffs for the various clusters
        it is not used if nl != None
        (plus one since the first cutoff is for pairs)
    """

    bothways = False
    lattice_neighbors = []
    # get structure
    if structure is None:
        if atoms is None:
            raise Exception(
                "Error: both atoms and structure is None in get_mbnls")
        else:
            structure = structure_from_atoms(atoms)

    # get neigbhorlists
    if neighborlists is None:
        neighborlists = []
        if cutoffs is None:
            raise Exception(
                "Error: both nl and cutoffs is None in count clusters")
        else:
            neighborlists = get_neighborlists(
                structure=structure, cutoffs=cutoffs)
    else:
        # build the neighborlists
        for nl in neighborlists:
            nl.build(structure)

    order = len(cutoffs) + 1

    mbnl = ManybodyNeighborlist()
    # add the mbnl lattice neighbors
    if order >= 2:
        for lattice_index in range(structure.size()):
            lattice_neighbor = mbnl.build(neighborlists, lattice_index, bothways)
            #print("order: {} len of  mbnl {}".format(order, len(lattice_neighbor)))
            for lat_nbrs in lattice_neighbor:
                lattice_neighbors.append(lat_nbrs)

    # add the pairs and singlets
    for lattice_index in range(structure.size()):
        lat_nbr_i = LatticeNeighbor(lattice_index, [.0, .0, .0])
        lattice_neighbors.append(([lat_nbr_i], []))  # singlet
        lattice_neighbors.append(
            ([lat_nbr_i], neighborlists[0].get_neighbors(lattice_index)))

    return lattice_neighbors
