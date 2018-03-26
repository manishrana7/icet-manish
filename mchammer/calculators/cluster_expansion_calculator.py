from mchammer.calculators.base_calculator import BaseCalculator


class ClusterExpansionCalculator(BaseCalculator):
    """
    Cluster expansion calculator.

    Class for efficiently calculating the
    cluster expansion property
    for a specific structure

    Parameters
    ----------
    atoms : ASE Atoms object
        the structure which the calculator can use
        to optimize the calculate function.

    cluster_expansion : icet ClusterExpansion object
    name : str
        Human readable identifier for this calculator

    """

    def __init__(self, atoms, cluster_expansion,
                 name='Cluster Expansion Calculator'):
        super().__init__(atoms=atoms, name=name)

        self._cluster_expansion = cluster_expansion

    @property
    def cluster_expansion(self):
        """
        icet cluster expansion.
        """
        return self._cluster_expansion

    def calculate_total(self):
        """
        Calculator the total property of the current configuration.

        Return:
            total_property : float
        """
        return self.cluster_expansion.predict(self.atoms)

    def calculate_local_contribution(self, indices):
        """
        return the local contribution of the index in indices.

        Return:
            local_contribution : float
        """
        local_contribution = 0
        for index in indices:
            local_contribution += self._calculate_local_contribution(index)

        return local_contribution

    def _calculate_local_contribution(self, index):
        """
        Internal method to calculate the local contribtuion for one
        index.
        """
        return self.calculate_total()
