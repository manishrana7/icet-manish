from mchammer.observers.base_observer import BaseObserver


class ClusterExpansionObserver(BaseObserver):
    """
    Cluster expansion observer.

    Parameters
    ----------
    cluster_expansion : :class:`icet:ClusterExpansion`
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

    def __init__(self, cluster_expansion, interval=None,
                 tag='ClusterExpansionObserver'):
        super().__init__(interval=interval, return_type=float, tag=tag)
        self._cluster_expansion = cluster_expansion

        if interval is None:
            raise ValueError("The value of interval must be specified")

    def get_observable(self, atoms) -> float:
        """
        Returns the value of the property from a cluster expansion model
        for a given atomic configuration.

        Parameters
        ----------
        atoms : :class:`ase:Atoms`
            input atomic structure.
        """
        return self._cluster_expansion.predict(atoms) * len(atoms)
