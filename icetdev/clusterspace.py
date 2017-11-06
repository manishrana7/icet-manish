from _icetdev import ClusterSpace
from icetdev import Structure
from icetdev.orbit_list import create_orbit_list
from ase import Atoms


def create_clusterspace(atoms, cutoffs, chemical_symbols,
                        Mi=None, verbosity=0):
    '''
    Creates a clusterspace.

    Parameters
    ----------
    atoms : ASE atoms object / icet structure object (bi-optional)
        The structure on which to base the clusterspace on.
    cutoffs : list of floats
        cutoff radii per order that define the clusterspace
    chemical_symbols : list of strings
        list of chemical symbols, each of which must map to an element of the
        periodic table
    verbosity : int
        verbosity level

    Todo
    ----
    * document `Mi`
    '''

    if isinstance(atoms, Atoms):
        structure = Structure.from_atoms(atoms)
    elif isinstance(atoms, Structure):
        structure = atoms
    else:
        msg = 'Unknown structure format'
        msg += ' {} (create_clusterspace)'.format(type(atoms))
        raise Exception(msg)

    orbit_list = create_orbit_list(structure, cutoffs, verbosity=verbosity)
    orbit_list.sort()
    if Mi is None:
        Mi = len(chemical_symbols)
    if isinstance(Mi, dict):
        Mi = get_Mi_from_dict(Mi, orbit_list.get_primitive_structure())
    if not isinstance(Mi, list):
        if not isinstance(Mi, int):
            raise Exception('Mi has wrong type in create_clusterspace')
        else:
            Mi = [Mi] * len(orbit_list.get_primitive_structure())

    msg = ['len(Mi) does not equal len(primitive_structure).']
    msg += ['{} != {}'.format(len(Mi),
                              len(orbit_list.get_primitive_structure()))]
    assert len(Mi) == len(orbit_list.get_primitive_structure()), ' '.join(msg)
    clusterspace = ClusterSpace(Mi, chemical_symbols, orbit_list)
    clusterspace.cutoffs = cutoffs
    return clusterspace


def __clusterspace__repr__(self):
    '''
    String representation of the clusterspcace
    '''
    s = ['Clusterspace']
    s += ['Subelements: {}'.format(' '.join(self.get_chemical_symbols()))]
    s += ['Cutoffs: {}'.format(' '.join(['{}'.format(co)
                                         for co in self.cutoffs]))]
    s += ['Total size of clusterspace: {}'.format(len(self))]

    t = ['Cluster bodies']
    t += ['cluster radius']
    t += ['multiplicity']
    t += ['(clusterspace index']
    t += [':orbit index)']
    t += ['mc vector']
    s += [' : '.join(t)]
    s += ['{:-^{l}}'.format(l=len(' : '.join(t)))]

    i = 0
    print_threshold = 50
    while i < len(self):
        if len(self) > print_threshold and i > 10 and i < len(self) - 10:
            i = len(self) - 10
            s += ['......\n\n']

        clusterspace_info = self.get_clusterspace_info(i)
        orbit_index = clusterspace_info[0]
        mc_vector = clusterspace_info[1]

        cluster = self.get_orbit(orbit_index).get_representative_cluster()
        multiplicity = len(self.get_orbit(orbit_index).get_equivalent_sites())

        t = ['{}'.format(len(cluster))]
        t += ['{.4f}'.format(cluster.get_geometrical_size())]
        t += ['{}'.format(multiplicity)]
        t += ['({}'.format(i)]
        t += ['{})'.format(orbit_index)]
        t += ['{}'.format(mc_vector)]
        s += [' : '.join(t)]

        i += 1

    return '\n'.join(s)


ClusterSpace.__repr__ = __clusterspace__repr__


def get_singlet_info(crystal_structure, return_clusterspace=False):
    '''
    Get information about the singlets in this structure.

    crystal_structure: icet structure object or ASE atoms object

    return_clusterspace: bool
        If true it will return the created clusterspace
    '''

    # create dummy elements and cutoffs
    subelements = ['H', 'He']
    cutoffs = [0.0]

    clusterspace = create_clusterspace(crystal_structure, cutoffs, subelements)

    singlet_data = []

    for i in range(len(clusterspace)):
        clusterspace_info = clusterspace.get_clusterspace_info(i)
        orbit_index = clusterspace_info[0]
        cluster = clusterspace.get_orbit(
            orbit_index).get_representative_cluster()
        multiplicity = len(clusterspace.get_orbit(
            orbit_index).get_equivalent_sites())
        assert len(cluster) == 1, \
            'Cluster in singlet only clusterspace has non-singlets'

        singlet = {}
        singlet['orbit_index'] = orbit_index
        singlet['sites'] = clusterspace.get_orbit(
            orbit_index).get_equivalent_sites()
        singlet['multiplicity'] = multiplicity
        singlet['representative_site'] = clusterspace.get_orbit(
            orbit_index).get_representative_sites()
        singlet_data.append(singlet)

    if return_clusterspace:
        return singlet_data, clusterspace
    else:
        return singlet_data


def view_singlets(structure):
    '''
    Visualize the singlets in the structure,
    singlet 0 is represented by a H,
    singlet 1 is represented by a He etc...
    '''

    cluster_data, clusterspace = get_singlet_info(structure,
                                                  return_clusterspace=True)

    primitive_atoms = clusterspace.get_primitive_structure().to_atoms()

    from ase.visualize import view
    from ase.data import chemical_symbols

    for singlet in cluster_data:
        for site in singlet['sites']:
            element = chemical_symbols[singlet['orbit_index']]
            atom_index = site[0].index
            primitive_atoms[atom_index].symbol = element

    view(primitive_atoms)


def get_Mi_from_dict(Mi, structure):
    '''
    Mi maps orbit index to allowed components
    this function will return a list, where
    Mi_ret[i] will be the allowed components on site index i

    '''
    cluster_data = get_singlet_info(structure)
    Mi_ret = [-1] * len(structure)
    for singlet in cluster_data:
        for site in singlet['sites']:
            Mi_ret[site[0].index] = Mi[singlet['orbit_index']]

    for all_comp in Mi_ret:
        if all_comp == -1:
            raise Exception(
                'Error: the calculated Mi from dict did not cover all sites on input structure. \n Were all sites in primitive mapped?')

    return Mi_ret


def __get_clustervector(self, crystal_structure):
    '''
    Get clustervector

    crystal_structure: either an ASE Atoms object or icet structure


    '''

    if isinstance(crystal_structure, Atoms):
            structure = Structure.from_atoms(crystal_structure)
    elif not isinstance(crystal_structure, Structure):
        print(type(crystal_structure))
        raise Exception(
            'Error: no known crystal structure format added to function create_clusterspace')
    else:
        structure = crystal_structure
    return self._get_clustervector(structure)

ClusterSpace.get_clustervector = __get_clustervector
