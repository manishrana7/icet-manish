from _icetdev import NeighborList
from .structure import Structure
from ase import Atoms


def get_neighbor_lists(atoms, cutoffs=None):
    '''
    Return list of icet neigbhorlists from a configuration and cutoffs

    Parameters
    ----------
    atoms : ASE Atomss object / icet Structure object (bi-optional)
        atomic configuration
    cutoffs:
        positive floats indicating the cutoffs for the various clusters

    Returns
    -------
    list of neighbor_lists
    '''

    # deal with different types of structure objects
    if isinstance(atoms, Atoms):
        structure = Structure.from_atoms(atoms)
    elif isinstance(atoms, Structure):
        structure = atoms
    else:
        msg = ['Unknown structure format']
        msg += ['{} (ClusterSpace)'.format(type(atoms))]
        raise Exception(' '.join(msg))

    neighbor_lists = []
    if cutoffs is None:
        raise Exception('Both n and cutoffs is None in count clusters')
    else:
        for co in cutoffs:
            nl = NeighborList(co)
            neighbor_lists.append(nl)

    # build the neighbor_lists
    for nl in neighbor_lists:
        nl.build(structure)

    return neighbor_lists
