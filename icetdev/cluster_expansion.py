import numpy as np
from icetdev import clusterspace

class ClusterExpansion(object):
    """
        Cluster expansion model

        attributes:
        -----------
        clusterspace : icet clusterspace object
        parameters : list of floats

    """

    def __init__(self, clusterspace, parameters):
        """
        Init
        Parameters:
        -----------
            :param clusterspace: icet ClusterSpace object
            :param parameters: list of floats
        """
        self._clusterspace = clusterspace
        self._parameters = parameters

    def predict(self, structure):
        """
        Predict a property from the input structure

        Parameters:
        -----------
            :param structure: icet Structure or ASE Atoms object

        Returns:
        -------
        float
            The predicted property of the structure
        """

        clustervector = self.clusterspace.get_clustervector(structure)        
        property = np.dot(clustervector, self.parameters)
        return property

    @property
    def clusterspace(self):
        '''icet clusterspace object'''
        return self._clusterspace

    @property
    def parameters(self):
        '''parameters'''
        return self._parameters
