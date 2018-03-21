from abc import ABC, abstractmethod
from mchammer.data_container import DataContainer


class BaseEnsemble(ABC):
    """
    Base ensemble abstract class.

    Attributes
    ----------
    name : str
        Name of the ensemble.
    accepted_trials : int
        Number of accepted trial steps.
    total_trials : int
        Numer of total trial steps.
    """

    def __init__(self, atoms, calculator, name='BaseEnsemble',
                 data_container=None, random_seed=42):

        self._calculator = calculator
        self._atoms = atoms
        self._name = name

        self.accepted_trials = 0
        self.total_trials = 0
        self._random_seed = random_seed
        self._boltzmann_constant = 8.6173303e-5
        self._observers = {}

        if data_container is None:
            self._data_container = DataContainer()
        else:
            pass

    @property
    def structure(self):
        """
        ASE Atoms object :
        The current state of the  structure being sampled in the ensemble.
        """
        return self._atoms

    @property
    def data_container(self):
        """
        mchammer DataContainer.
        """
        return self._data_container

    @property
    def observers(self):
        """
        dict of mchammer Observer's.
        """
        return self._observers

    @property
    def acceptance_ratio(self):
        """
        Returns the acceptance ratio,
        i.e. accepted_trials / total_trials.
        """
        return self.accepted_trials / self.total_trials

    @property
    def calculator(self):
        """
        mchammer Calculator.
        """
        return self._calculator

    def run(self, trial_moves, observation_interval=None):
        """
        Sample the ensemble `trial_moves` steps.

        Parameters:
        -----------

        trial_moves: int
            number of steps to run in total
        """
        pass

    @abstractmethod
    def do_trial_move(self):
        pass

    @property
    def name(self):
        """
        str : Name of BaseEnsemble.
        """
        return self._name

    @property
    def random_seed(self):
        """
        int : Random seed used in random number generator.
        """
        return self._random_seed

    def next_random_number(self):
        """
        Return the next random number from the RNG.
        """
        pass

    @property
    def kB(self):
        """
        Boltzmann constant in eV / K.
        """
        return self._boltzmann_constant

    @abstractmethod
    def _run(self, trial_moves):
        """
        Internal running of the MC without interruption.

        trial_moves : int
           Number of trial moves to run without stopping.
        """
        pass

    def attach_observer(self, observer, tag=None):
        """
        Attaches an observer to the ensemble.

        Parameters:
        ----------
        observer : mchammer Observer object
        """

        if tag is not None:
            self.observers[tag] = observer
        else:
            self.observers[observer.tag] = observer
