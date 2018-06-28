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
        cluster expansion model to be used for observation
    tag : str (default: `ClusterExpansionObserver`)
        human readable observer name
    interval : int
        interval at which observation is done during the Monte Carlo
        simulation

    Attributes
    ----------
    tag : str
        name of observer
    interval : int
        observation interval
    """

    def __init__(self, atoms, cluster_expansion, scaling=None,
                 interval=None, tag='ClusterExpansionObserver'):
        super().__init__(interval=interval, return_type=float, tag=tag)
        self._cluster_expansion = cluster_expansion
        self.atoms = atoms

        if interval is None:
            raise ValueError("The value of interval must be specified")

        if scaling is None:
            self._property_scaling = len(atoms)
        else:
            self._property_scaling = scaling

    def get_observable(self, occupations: List[int]) -> float:
        """
        Returns the value of a property from a cluster expansion model.

        Parameters
        ----------
        occupations
            list of atomic species or vector occupation of the atomic
            configuration.
        """
        self.atoms.set_atomic_numbers(occupations)
        return self.cluster_expansion.predict(self.atoms) * \
            self._property_scaling
