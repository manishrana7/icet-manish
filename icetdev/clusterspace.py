from _icetdev import ClusterSpace
from icetdev.structure import structure_from_atoms
from icetdev.orbitList import create_orbit_list


def create_clusterspace(subelements, cutoffs, atoms=None, structure=None,Mi=None, verbosity=0):
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

    orbitList = create_orbit_list(structure, cutoffs, verbosity=verbosity)
    orbitList.sort()
    if Mi == None:
        Mi = len(subelements)
    if isinstance(Mi,dict):
        Mi = get_Mi_from_dict(Mi, orbitList.get_primitive_structure())
    if not isinstance(Mi,list):        
        if not isinstance(Mi, int):
            raise Exception("Error: Mi has wrong type in create_clusterspace")
        else:
            Mi = [Mi] * len(orbitList.get_primitive_structure())
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
        rep += str(co) + " "
    rep += " \n"
    rep += " Total number of dimensions {} \n".format(len(self))
    rep += "clusterspace index : orbit index : Cluster order : cluster radius : multiplicity  : mc vector\n"
    rep += "-------------------------------------\n"
    i = 0
    print_threshold = 50
    while i < len(self):
        if len(self) > print_threshold and i > 10 and i < len(self) - 10:
            i = len(self) - 10
            rep += "\n...... \n\n"
        rep += "{}: ".format(i)
        clusterspace_info = self.get_clusterspace_info(i)
        orbit_index = clusterspace_info[0]
        mc_vector = clusterspace_info[1]

        cluster = self.get_orbit(orbit_index).get_representative_cluster()
        multiplicity = len(self.get_orbit(orbit_index).get_equivalent_sites())

        rep += " {0} : {1} {2:.4f} : {3} : {4} ".format(len(cluster),
                                                        orbit_index,
                                                        cluster.get_geometrical_size(), multiplicity, mc_vector)

        rep += "\n"
        i += 1
    return rep


ClusterSpace.__repr__ = __represent_clusterspace


def get_singlet_info(structure=None, atoms=None, return_clusterspace=False):
    """
    Get information about the singlets in this structure.

    structure: icet structure object
        either this or atoms needs to be defined
    atoms: ASE atoms object
        either this or structure needs to be defined
    return_clusterspace: bool
        If true it will return the created clusterspace        
    """

    # if structure==None:
    #     if atoms == None:
    #         Exception("Error: neither structure or atoms were defined in function: Clusterspace.py : get_singlet_info")
    #     else:
    #         structure = structure_from_atoms(atoms)

    # create dummy elements and cutoffs
    subelements = ["H", "He"]
    cutoffs = [0.0]

    clusterspace = create_clusterspace(
        subelements, cutoffs, atoms=atoms, structure=structure)

    singlet_data = []

    for i in range(len(clusterspace)):
        clusterspace_info = clusterspace.get_clusterspace_info(i)
        orbit_index = clusterspace_info[0]
        cluster = clusterspace.get_orbit(
            orbit_index).get_representative_cluster()
        multiplicity = len(clusterspace.get_orbit(
            orbit_index).get_equivalent_sites())
        assert len(
            cluster) == 1, "Error: cluster in singlet only clusterspace has non-singlets"

        singlet = {}
        singlet["orbit_index"] = orbit_index
        singlet["sites"] = clusterspace.get_orbit(
            orbit_index).get_equivalent_sites()
        singlet["multiplicity"] = multiplicity
        singlet["representative_site"] = clusterspace.get_orbit(
            orbit_index).get_representative_sites()
        singlet_data.append(singlet)

    if return_clusterspace:
        return singlet_data, clusterspace
    else:
        return singlet_data


def view_singlets(structure=None, atoms=None):
    """
    Visualize the singlets in the structure,
    singlet 0 is represented by a H,
    singlet 1 is represented by a He etc...
    """

    cluster_data, clusterspace = get_singlet_info(
        structure=structure, atoms=atoms, return_clusterspace=True)

    primitive_atoms = clusterspace.get_primitive_structure().to_atoms()

    from ase.visualize import view
    from ase.data import atomic_numbers, chemical_symbols

    for singlet in cluster_data:
        for site in singlet["sites"]:
            element = chemical_symbols[singlet["orbit_index"]]
            atom_index = site[0].index
            primitive_atoms[atom_index].symbol = element

    view(primitive_atoms)




def get_Mi_drom_dict(Mi, structure):
    cluster_data = get_singlet_info(atoms=structure)
    Mi_ret = [-1] * len(structure)
    for singlet in cluster_data:
        for site in singlet["sites"]:
            Mi_ret[ site[0].index ] = Mi[singlet["orbit_index"]]

    return Mi_ret