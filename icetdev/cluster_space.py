from _icetdev import ClusterSpace as _ClusterSpace
from icetdev import Structure
from icetdev.orbit_list import create_orbit_list
from icetdev.permutation_map import vacuum_on_non_pbc
from ase import Atoms
import numpy as np


class ClusterSpace(_ClusterSpace):

    def __init__(self, atoms, cutoffs, chemical_symbols,
                 Mi=None, verbosity=0):
        '''
        Prepare a cluster space.

        Parameters
        ----------
        atoms : ASE Atoms object / icet Structure object (bi-optional)
            atomic configuration
        cutoffs : list of floats
            cutoff radii per order that define the cluster space
        chemical_symbols : list of strings
            list of chemical symbols, each of which must map to an element of
            the periodic table
        verbosity : int
            verbosity level

        Todo
        ----
        * document `Mi`
        '''

        # deal with different types of structure objects
        if isinstance(atoms, Atoms):
            structure = Structure.from_atoms(atoms)
        elif isinstance(atoms, Structure):
            structure = atoms
        else:
            msg = 'Unknown structure format'
            msg += ' {} (ClusterSpace)'.format(type(atoms))
            raise Exception(msg)

        # set up orbit list
        orbit_list = create_orbit_list(structure, cutoffs, verbosity=verbosity)
        orbit_list.sort()
        if Mi is None:
            Mi = len(chemical_symbols)
        if isinstance(Mi, dict):
            Mi = get_Mi_from_dict(Mi, orbit_list.get_primitive_structure())
        if not isinstance(Mi, list):
            if not isinstance(Mi, int):
                raise Exception('Mi has wrong type (ClusterSpace)')
            else:
                Mi = [Mi] * len(orbit_list.get_primitive_structure())
        msg = ['len(Mi) does not equal len(primitive_structure);']
        msg += ['{} != {}'.format(len(Mi),
                                  len(orbit_list.get_primitive_structure()))]
        msg = ' '.join(msg)
        assert len(Mi) == len(orbit_list.get_primitive_structure()), msg

        # call (base) C++ constructor
        _ClusterSpace.__init__(self, Mi, chemical_symbols, orbit_list)
        self.cutoffs = cutoffs

    def get_overview(self, print_threshold=None, parameters=None):
        '''
        Obtain an overview of the cluster space in terms of the orbits (order,
        radius, multiplicity etc) including optionally the associated
        parameters (i.e., effective cluster interactions=ECIs).

        Parameters
        ----------
        print_threshold : int
            if the number of orbits exceeds this number print dots
        parameters : list
            parameters (ECIs) to be printed along side orbit information

        Returns
        -------
        multi-line string
            string representation of the cluster space.
        '''

        def repr_cluster(index, cluster, multiplicity=0,
                         orbit_index=0, mc_vector=[0] * 5, param=None):
            from collections import OrderedDict
            fields = OrderedDict([
                ('order',   '{:2}'.format(cluster.order)),
                ('radius',  '{:9.4f}'.format(cluster.geometrical_size)),
                ('multiplicity', '{:4}'.format(multiplicity)),
                ('index',     '{:4}'.format(index)),
                ('orbit',     '{:4}'.format(orbit_index))])
            if param is not None:
                fields['ECI'] = '{:13.6e}'.format(param)
            fields['MC vector'] = '{:}'.format(mc_vector)
            s = []
            for name, value in fields.items():
                n = max(len(name), len(value))
                if index < 0:
                    s += ['{s:^{n}}'.format(s=name, n=n)]
                else:
                    s += ['{s:^{n}}'.format(s=value, n=n)]
            return ' | '.join(s)

        if parameters is not None:
            if len(parameters) - 1 != len(self):
                raise Exception('Length of parameter vector incompatible'
                                ' with cluster space size')
            param = 1.0
        else:
            param = None

        # basic information
        cluster = self.get_orbit(0).get_representative_cluster()
        n = len(repr_cluster(-1, cluster, param=param))
        s = []
        s += ['{s:-^{n}}'.format(s=' Cluster Space ', n=n)]
        s += [' subelements: {}'.format(' '.join(self.get_atomic_numbers()))]
        s += [' cutoffs: {}'.format(' '.join(['{}'.format(co)
                                              for co in self.cutoffs]))]
        s += [' number of orbits: {}'.format(len(self))]

        # table header
        horizontal_line = '{s:-^{n}}'.format(s='', n=n)
        s += [horizontal_line]
        s += [repr_cluster(-1, cluster, param=param)]
        s += [horizontal_line]

        # table body
        index = 0
        while index < len(self):
            if (print_threshold is not None and
                    len(self) > print_threshold and
                    index > 10 and
                    index < len(self) - 10):
                index = len(self) - 10
                s += [' ...']

            cluster_space_info = self.get_cluster_space_info(index)
            orbit_index = cluster_space_info[0]
            mc_vector = cluster_space_info[1]
            cluster = self.get_orbit(orbit_index).get_representative_cluster()
            multiplicity = len(self.get_orbit(
                               orbit_index).get_equivalent_sites())
            if parameters is not None:
                param = parameters[index+1]
            s += [repr_cluster(index, cluster, multiplicity,
                               orbit_index, mc_vector, param)]

            index += 1

        s += [horizontal_line]

        return '\n'.join(s)

    def __repr__(self):
        '''
        String representation of the cluster space.
        '''
        return self.get_overview(50)

    def get_cluster_vector(self, atoms):
        '''
        Returns the cluster vector for a structure.

        Parameters
        ----------
        atoms : ASE Atoms object / icet Structure object (bi-optional)
            atomic configuration

        Returns
        -------
        NumPy array
            the cluster vector
        '''
        if isinstance(atoms, Atoms):
            structure = Structure.from_atoms(atoms)
        elif isinstance(atoms, Structure):
            structure = atoms
        else:
            msg = 'Unknown structure format'
            msg += ' {} (ClusterSpace.get_cluster_vector)'.format(type(atoms))
            raise Exception(msg)

        # if pbc is not true one needs to massage the structure a bit
        if not np.array(structure.get_pbc()).all():
            atoms = structure.to_atoms()
            vacuum_on_non_pbc(atoms)
            structure = Structure.from_atoms(atoms)
        else:
            atoms = structure.to_atoms()
            atoms.wrap()
            structure = Structure.from_atoms(atoms)
        return _ClusterSpace.get_cluster_vector(self, structure)


def get_singlet_info(atoms, return_cluster_space=False):
    '''
    Retrieve information concerning the singlets in the input structure.

    Parameters
    ----------
    atoms : ASE Stoms object / icet Structure object (bi-optional)
        atomic configuration
    return_cluster_space : boolean
        return the cluster space created during the process

    Returns
    -------
    list of dicts
        each dictionary in the list represents one orbit
    ClusterSpace object (optional)
        cluster space created during the process
    '''

    # create dummy elements and cutoffs
    subelements = ['H', 'He']
    cutoffs = [0.0]

    cs = ClusterSpace(atoms, cutoffs, subelements)

    singlet_data = []

    for i in range(len(cs)):
        cluster_space_info = cs.get_cluster_space_info(i)
        orbit_index = cluster_space_info[0]
        cluster = cs.get_orbit(orbit_index).get_representative_cluster()
        multiplicity = len(cs.get_orbit(orbit_index).get_equivalent_sites())
        assert len(cluster) == 1, 'Cluster space contains not only singlets'

        singlet = {}
        singlet['orbit_index'] = orbit_index
        singlet['sites'] = cs.get_orbit(orbit_index).get_equivalent_sites()
        singlet['multiplicity'] = multiplicity
        singlet['representative_site'] = cs.get_orbit(
            orbit_index).get_representative_sites()
        singlet_data.append(singlet)

    if return_cluster_space:
        return singlet_data, cs
    else:
        return singlet_data


def get_singlet_configuration(atoms, to_primitive=False):
    '''
    Return atomic configuration decorated with a different element for each
    Wyckoff site. This is useful for visualization and analysis.

    Parameters
    ----------
    atoms : ASE Atoms object / icet Structure object (bi-optional)
        atomic configuration
    to_primitive : boolean
        if True the input structure will be reduced to its primitive unit cell
        before processing

    Returns
    -------
    ASE Atoms object
        structure with singlets highlighted by different elements
    '''

    from ase.data import chemical_symbols
    cluster_data, cluster_space = get_singlet_info(atoms,
                                                   return_cluster_space=True)

    if to_primitive:
        singlet_configuration = \
            cluster_space.get_primitive_structure().to_atoms()
        for singlet in cluster_data:
            for site in singlet['sites']:
                element = chemical_symbols[singlet['orbit_index'] + 1]
                atom_index = site[0].index
                singlet_configuration[atom_index].symbol = element
    else:
        singlet_configuration = atoms.copy()
        singlet_configuration = vacuum_on_non_pbc(singlet_configuration)
        orbit_list = cluster_space.get_orbit_list()
        orbit_list_supercell \
            = orbit_list.get_supercell_orbit_list(singlet_configuration)
        for singlet in cluster_data:
            for site in singlet['sites']:
                element = chemical_symbols[singlet['orbit_index'] + 1]
                sites = orbit_list_supercell.get_orbit(
                    singlet['orbit_index']).get_equivalent_sites()
                for lattice_site in sites:
                    k = lattice_site[0].index
                    singlet_configuration[k].symbol = element

    return singlet_configuration


def view_singlets(atoms, to_primitive=False):
    '''
    Visualize singlets in a structure using the ASE graphical user interface.

    Parameters
    ----------
    atoms : ASE Atoms object / icet Structure object (bi-optional)
        atomic configuration
    to_primitive : boolean
        if True the input structure will be reduced to its primitive unit cell
        before processing
    '''
    from ase.visualize import view
    singlet_configuration = get_singlet_configuration(
        atoms, to_primitive=to_primitive)
    view(singlet_configuration)


def get_Mi_from_dict(Mi, structure):
    '''
    Mi maps orbit index to allowed components
    this function will return a list, where
    Mi_ret[i] will be the allowed components on site index i

    Parameters
    ----------
    Mi : xx
        xx
    atoms : ASE Atoms object / icet Structure object (bi-optional)
        atomic configuration

    Returns
    -------
    list
        xxx

    Todo
    ----
    * rename function, remove `Mi`
    * complete docstring
    '''
    cluster_data = get_singlet_info(structure)
    Mi_ret = [-1] * len(structure)
    for singlet in cluster_data:
        for site in singlet['sites']:
            Mi_ret[site[0].index] = Mi[singlet['orbit_index']]

    for all_comp in Mi_ret:
        if all_comp == -1:
            s = ['The calculated Mi from dict did not cover all sites of'
                 ' the input structure.']
            s += ['Were all sites in primitive mapped?']
            raise Exception('\n'.join(s))

    return Mi_ret
