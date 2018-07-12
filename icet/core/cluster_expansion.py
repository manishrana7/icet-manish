import numpy as np
import pickle
import re
from icet import ClusterSpace


class ClusterExpansion(object):
    '''
    Cluster expansion model

    Attributes
    ----------
    cluster_space : ClusterSpace object
        cluster space that was used for constructing the cluster expansion
    parameters : list of floats
        effective cluster interactions (ECIs)
    '''

    def __init__(self, cluster_space, parameters):
        '''
        Initialize a ClusterExpansion object.

        Parameters
        ----------
        cluster_space : ClusterSpace object
            the cluster space to be used for constructing the cluster expansion
        parameters : list of floats
            effective cluster interactions (ECIs)
        '''
        if len(cluster_space) != len(parameters):
            raise ValueError('cluster_space and parameters must have the same'
                             ' length ({} != {})'.format(len(cluster_space),
                                                         len(parameters)))
        self._cluster_space = cluster_space
        self._parameters = parameters

    def predict(self, structure):
        '''
        Predict the property of interest (e.g., the energy) for the input
        structure using the cluster expansion.

        Parameters
        ----------
        structure : ASE Atoms object / icet Structure (bi-optional)
            atomic configuration

        Returns
        -------
        float
            property value predicted by the cluster expansion
        '''
        cluster_vector = self.cluster_space.get_cluster_vector(structure)
        prop = np.dot(cluster_vector, self.parameters)
        return prop

    @property
    def cluster_space(self):
        '''ClusterSpace object : cluster space the cluster expansion is
        based on'''
        return self._cluster_space

    @property
    def parameters(self):
        '''list of floats : effective cluster interactions (ECIs)'''
        return self._parameters

    def write(self, filename):
        """
        Write Cluster expansion to file.

        Parameters
        ---------
        filename : str with filename to saved
        cluster space.
        """

        self.cluster_space.write(filename)

        with open(filename, 'rb') as handle:
            data = pickle.load(handle)

        data['parameters'] = self.parameters

        with open(filename, "wb") as handle:
            pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def read(filename):
        """
        Read cluster expansion from file.

        Parameters
        ---------
        filename : str with filename to saved
        cluster space.
        """
        cs = ClusterSpace.read(filename)
        if isinstance(filename, str):
            with open(filename, 'rb') as handle:
                data = pickle.load(handle)
        else:
            with open(filename) as handle:
                data = pickle.load(handle)
        parameters = data['parameters']

        return ClusterExpansion(cs, parameters)

    def __len__(self) -> int:
        return len(self._parameters)

    def _get_string_representation(self, print_threshold=None,
                                   print_minimum=10):
        """
        String representation of the cluster expansion.
        """
        cluster_space_repr = self._cluster_space._get_string_representation(
            print_threshold, print_minimum).split('\n')
        # rescale width
        eci_col_width = max(
            len('{:4.1f}'.format(max(self._parameters, key=abs))), len('ECI'))
        width = len(cluster_space_repr[0]) + len(' | ') + eci_col_width

        s = []
        s += ['{s:=^{n}}'.format(s=' Cluster Expansion ', n=width)]
        s += [t for t in cluster_space_repr if re.search(':', t)]

        # table header
        s += [''.center(width, '-')]
        t = [t for t in cluster_space_repr if 'index' in t]
        t += ['{s:^{n}}'.format(s='ECI', n=eci_col_width)]
        s += [' | '.join(t)]
        s += [''.center(width, '-')]

        # table body
        index = 0
        while index < len(self):
            if (print_threshold is not None and
                    len(self) > print_threshold and
                    index >= print_minimum and
                    index <= len(self) - print_minimum):
                index = len(self) - print_minimum
                s += [' ...']
            pattern = r'^{:4}'.format(index)
            t = [t for t in cluster_space_repr if re.match(pattern, t)]
            t += ['{s:^{n}}'.format(s=self._parameters[index],
                                    n=eci_col_width)]
            s += [' | '.join(t)]
            index += 1
        s += [''.center(width, '=')]

        return '\n'.join(s)

    def __repr__(self) -> str:
        """String representation."""
        return self._get_string_representation(print_threshold=50)

    def print_overview(self, print_threshold=None, print_minimum=10):
        """
        Print an overview of the cluster space in terms of the orbits (order,
        radius, multiplicity etc).

        Parameters
        ----------
        print_threshold : int
            if the number of orbits exceeds this number print dots
        print_minimum : int
            number of lines printed from the top and the bottom of the orbit
            list if `print_threshold` is exceeded
        """
        print(self._get_string_representation(print_threshold=print_threshold,
                                              print_minimum=print_minimum))
