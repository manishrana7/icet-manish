from mchammer.observers.base_observer import BaseObserver
from typing import List


class ClusterExpansionObserver(BaseObserver):
    """
    Cluster expansion observer.

    Description of this class.

    Parameters
    ----------
    atoms : :class:`ase:Atoms`
        atomic structure to be used as a reference by the observer.
    cluster_expansion :class:`icet:ClusterExpansion` 
        
    tag : str (default: `ClusterExpansionObserver`)
        human readable observer name

    Attributes
    ----------
    tag : str
        name of observer
    interval : int
        observation interval
    """

    def __init__(self, atoms, cluster_expansion,
                 tag='ClusterExpansionObserver', interva),
                 atoms, cluster_expansion):
        super().__init__(interval=interval, return_type=return_type, tag=tag)
        self._cluster_expansion = cluster_expansion
        self.atom = atom

    def get_observable(self, occupations: List[int]) -> float:
        """
        Returns an observation.
        """
        self.atoms.set_atomic_numbers(occupations)
        return self.cluster_expansion.predict(self.atoms) * \
            self._property_scaling
