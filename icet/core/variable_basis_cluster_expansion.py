from icet import (ClusterExpansion,
                  VariableBasisClusterSpace)
import numpy as np
import pickle
import tempfile
import tarfile


class VariableBasisClusterExpansion(ClusterExpansion):
    cluster_space_type = VariableBasisClusterSpace

    def __init__(self, cluster_space: VariableBasisClusterSpace, parameters: np.array,
                 metadata: dict = None):
        """
        Initializes a VariableBasisClusterExpansion object.

        Parameters
        ----------
        cluster_space
            cluster space to be used for constructing the cluster expansion
        parameters
            parameter vector
        metadata : dict
            metadata dictionary, user-defined metadata to be stored together
            with cluster expansion. Will be pickled when CE is written to file.
            By default contains icet version, username, hostname and date.

        Raises
        ------
        ValueError
            if the length of parameters list does not match the total number
            of parameters expected by the cluster space
        """
        if not isinstance(cluster_space, VariableBasisClusterSpace):
            raise ValueError('VariableBasisClusterExpansion can only be initlaized '
                             'with a VariableBasisClusterSpace')
        if cluster_space.total_number_of_parameters != len(parameters):
            raise ValueError('cluster_space ({}) and parameters ({}) must have'
                             ' the same length'.format(len(cluster_space), len(parameters)))
        self._cluster_space = cluster_space.copy()
        if isinstance(parameters, list):
            parameters = np.array(parameters)
        self._parameters = parameters

        # add metadata
        if metadata is None:
            metadata = dict()
        self._metadata = metadata
        self._add_default_metadata()

