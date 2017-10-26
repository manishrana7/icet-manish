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
    clusterspace.cutoffs = cutoffs
    return clusterspace


def __size_of_clusterspace(self):
    """"
    Returns the size of the clusterspace, i.e. the length of a clustervector in this space
    """
    return self.get_clusterspace_size()


ClusterSpace.__len__ = __size_of_clusterspace


def __represent_clusterspace(self):
    """
    String representation of the clusterspcace
    """
    rep = "Clusterspace \n"

    rep += "Subelements: "
    for el in self.get_elements():
        rep += el
    rep += " \n"

    rep += "Cutoffs: "
    for co in self.cutoffs:
        rep += str(co) +" "
    rep += " \n"
    rep += " Total number of dimensions {} \n".format(len(self))
    rep += "clusterspace index : Cluster order : cluster radius : multiplicity  : mc vector\n" 
    rep += "-------------------------------------\n"
    for i in range(len(self)):
        rep += "{}: ".format(i)
        clusterspace_info = self.get_clusterspace_info(i)
        orbit_index = clusterspace_info[0]
        mc_vector = clusterspace_info[1]

        cluster = self.get_orbit(orbit_index).get_representative_cluster()
        multiplicity = len(self.get_orbit(orbit_index).get_equivalent_sites())

        rep += " {0} : {1:.4f} : {2} : {3} ".format(len(cluster),
        cluster.get_geometrical_size(), multiplicity, mc_vector)
                                
        rep += "\n"
    return rep
ClusterSpace.__repr__ = __represent_clusterspace