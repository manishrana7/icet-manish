from _icetdev import ClusterSpace
from icetdev.structure import structure_from_atoms
from icetdev.orbitList import create_orbit_list


def create_clusterspace(subelements, cutoffs, atoms=None, structure=None, verbosity=0):
    """
    Creates a clusterspace.

    subelements: list of strings
        The strings are required to be one of the short descriptions in the periodic table
    cutoffs: list of floats
        Cutoffs that define the clusterspace
    atoms: ASE atoms object (bi-optional)
        The structure on which to base the clusterspace on.
        Atleast one of structure or atoms need to be non-None
    structure: icet structure object (bi-optional)         
        The structure on which to base the clusterspace on.
        Atleast one of structure or atoms need to be non-None
    verbosity: int
        sets the verbosity level

    """

    # get structure
    if structure is None:
        if atoms is None:
            raise Exception(
                "Error: both atoms and structure is None in get_mbnls")
        else:
            structure = structure_from_atoms(atoms)

    Mi = len(subelements)

    orbitList = create_orbit_list(structure, cutoffs, verbosity=verbosity)
    orbitList.sort()
    clusterspace = ClusterSpace(Mi, subelements, orbitList)
    return clusterspace