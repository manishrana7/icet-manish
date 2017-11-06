from _icetdev import ManybodyNeighborlist
from icetdev.lattice_site import LatticeSite
from icetdev.neighborlist import get_neighborlists
from icetdev.structure import Structure


def get_all_lattice_neighbors(atoms=None, structure=None,
                              neighborlists=None, cutoffs=None):
    '''
    Generate lattice neighbors for a configuration.

    Parameters
    ----------
    atoms : ASE Atoms object
        input structure (bi-optional)
    structure :  icet Structure object
        input structure (bi-optional)
    neighborlists : array/list of icet NeighborList objects
        neighbor lists; if None neighbor lists will be created, which
        requires cutoffs to be specified
    cutoffs : list of floats
        cutoffs for the clusters of increasing order; only used if
        `neighborlists` is None; note that the first value in the list refers
        to pairs

    Returns
    -------
    list of tuples
        the first part of each tuple identifies the lattice neighbor, the
        second part is a list of neighbors (by index)
    '''

    bothways = False
    lattice_neighbors = []
    # get structure
    if structure is None:
        if atoms is None:
            raise Exception('Need to specify either atoms or structure')
        else:
            structure = Structure.from_atoms(atoms)

    # get neigbhor lists
    if neighborlists is None:
        neighborlists = []
        if cutoffs is None:
            raise Exception('Need to specify either neighbor list or cutoffs')
        else:
            neighborlists = get_neighborlists(
                structure=structure, cutoffs=cutoffs)
    else:
        # build the neighbor lists
        for nl in neighborlists:
            nl.build(structure)

    order = len(cutoffs) + 1

    mbnl = ManybodyNeighborlist()
    # add the mbnl lattice neighbors
    if order >= 2:
        for lattice_index in range(len(structure)):
            lattice_neighbor = mbnl.build(neighborlists,
                                          lattice_index, bothways)
            for lat_nbrs in lattice_neighbor:
                lattice_neighbors.append(lat_nbrs)

    # add the pairs and singlets
    for lattice_index in range(len(structure)):
        lat_nbr_i = LatticeSite(lattice_index, [0.0, 0.0, 0.0])
        lattice_neighbors.append(([lat_nbr_i], []))  # singlet
        lattice_neighbors.append(
            ([lat_nbr_i], neighborlists[0].get_neighbors(lattice_index)))

    return lattice_neighbors
