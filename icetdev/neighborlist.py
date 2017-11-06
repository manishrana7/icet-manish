from _icetdev import Neighborlist
from icetdev.structure import Structure


def get_neighborlists(atoms=None, structure=None, cutoffs=None):
    '''
    Return list of icet neigbhorlists from a configuration and cutoffs

    Parameters
    ----------
    atoms:
        ASE atoms object (bi-optional) need atleast this or structure
    structure:
        icet structure object (bi-optional) need atleast this or atoms
    cutoffs:
        positive floats indicating the cutoffs for the various clusters

    Returns
    -------
    list of neighborlists
    '''

    if structure is None:
        if atoms is None:
            raise Exception(
                'Error: both atoms and structure is None in get_mbnls')
        else:
            structure = Structure.from_atoms(atoms)

    neighborlists = []
    if cutoffs is None:
        raise Exception('Both n and cutoffs is None in count clusters')
    else:
        for co in cutoffs:
            nl = Neighborlist(co)
            neighborlists.append(nl)

    # build the neighborlists
    for nl in neighborlists:
        nl.build(structure)

    return neighborlists
