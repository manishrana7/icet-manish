from mchammer.observers.base_observer import BaseObserver
from typing import List


class ClusterExpansionObserver(BaseObserver):
    """
    Cluster expansion observer.

    Description of this class.

    Parameters
    ----------
    atoms : ASE atoms object

    cluster_expansion : icet cluster expansion

    interval : int

    tag : str

    Attributes
    ----------
    cluster_expansion : icet cluster expansion

    tag : str
        human readable tag used for identifying the observer

    interval : int
        the observation interval
    """

    def __init__(self,  cluster_expansion, interval, tag, return_type):
        super().__init__(interval=interval, return_type=return_type,
                         tag=tag)
        self._cluster_expansion = cluster_expansion

    def get_observable(self, occupations: List[int]) -> float:
        pass
