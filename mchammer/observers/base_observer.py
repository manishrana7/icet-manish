"""
Definition of the base observer class.
"""
from abc import ABC, abstractmethod


class BaseObserver(ABC):
    """
    Base Observer class.

    Attributes
    ----------
    tag : str
        Human readable tag used for identifying the observer.
    """

    def __init__(self, interval, tag='BaseObserver'):
        self.tag = tag
        self.interval = interval

    @property
    @abstractmethod
    def return_type(self):
        """
        Data type of the observed data.
        """
        raise NotImplementedError

    @abstractmethod
    def get_observable(self):
        """
        Method used for extracting data.

        Returns:
            self.return_type()

        When implementing this method use the
        following names for the following types
        of data:
        ASE atoms object : `atoms`.
        list of chemical elements : `elements`.
        icet structure object : `structure`.
        icet cluster expansion : `cluster_expansion`.
        mchammer calculator : `calculator`.
        """
        raise NotImplementedError
