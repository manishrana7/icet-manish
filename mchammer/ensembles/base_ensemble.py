import random
from time import time
from abc import ABC, abstractmethod
from math import gcd

from mchammer.data_container import DataContainer
from mchammer.configuration_manager import ConfigurationManager
from mchammer.observers.base_observer import BaseObserver

import numpy as np


class BaseEnsemble(ABC):
    """
    Base ensemble abstract class.

    Attributes
    ----------
    name : str
        name of the ensemble.
    accepted_trials : int
        number of accepted trial steps.
    total_trials : int
        number of total trial steps.
    """

    def __init__(self, calculator, name='BaseEnsemble',
                 data_container=None, data_container_write_period=np.inf,
                 random_seed=None):

        self._calculator = calculator
        self._name = name
        self.data_container_write_period = data_container_write_period
        self.accepted_trials = 0
        self.total_trials = 0
        self._observers = {}
        self._step = 0

        if random_seed is None:
            self._random_seed = random.randint(0, 1e16)
        else:
            self._random_seed = random_seed
        random.seed(a=self.random_seed)

        if data_container is None:
            self._data_container = DataContainer(atoms=None,
                                                 ensemble_name=name,
                                                 random_seed=random_seed)
        else:
            raise NotImplementedError
        self._data_container_filename = None

        strict_constraints = self.calculator.occupation_constraints
        sublattices = [[i for i in range(len(self.calculator.atoms))]]
        self.configuration = ConfigurationManager(
            self._calculator.atoms.numbers, strict_constraints, sublattices)

    @property
    def structure(self):
        """
        ASE Atoms object :
        The current state of the  structure being sampled in the ensemble.
        """
        return self.calculator.atoms.copy()

    @property
    def step(self) -> int:
        """
        int : current MC trial step.
        """
        return self._step

    @property
    def data_container(self):
        """
        mchammer DataContainer object
        """
        return self._data_container

    @property
    def observers(self):
        """
        dict : mchammer observers.
        """
        return self._observers

    @property
    def acceptance_ratio(self):
        """
        float : the acceptance ratio,
        i.e. accepted_trials / total_trials.
        """
        return self.accepted_trials / self.total_trials

    @property
    def calculator(self):
        """
        mchammer Calculator.
        """
        return self._calculator

    def run(self, number_of_trial_steps: int, reset_step=False):
        """
        Sample the ensemble for `number_of_trial_steps` steps.

        Parameters:
        -----------
        number_of_trial_steps: int
            number of steps to run in total
        reset_step : bool
            if True the MC trial step counter will be initialized to zero
        """

        last_write_time = time()
        if reset_step:
            initial_step = 0
            final_step = number_of_trial_steps
            self._step = 0
        else:
            initial_step = self._step
            final_step = self._step + number_of_trial_steps
            # run mc so that we start at an interval which lands
            # on the observers interval
            if not initial_step == 0:
                first_run_interval = self.minimum_observation_interval -\
                    (initial_step -
                     (initial_step // self.minimum_observation_interval) *
                        self.minimum_observation_interval)
                first_run_interval = min(
                    first_run_interval, number_of_trial_steps)
                self._run(first_run_interval)
                initial_step += first_run_interval
                self._step += first_run_interval

        step = initial_step
        while step < final_step:
            uninterrupted_steps = min(
                self.minimum_observation_interval, final_step - step)
            if self._step % self.minimum_observation_interval == 0:
                self._observe_configuration(self._step)
            if self._data_container_filename is not None and \
                    time() - last_write_time > \
                    self.data_container_write_period:
                self.data_container.write(self._data_container_filename)

            self._run(uninterrupted_steps)
            step += uninterrupted_steps
            self._step += uninterrupted_steps

        # If we end on an observation interval we also observe
        if self._step % self.minimum_observation_interval == 0:
            self._observe_configuration(self._step)

        if self._data_container_filename is not None:
            self.data_container.write(self._data_container_filename)

    def _run(self, number_of_trial_steps):
        """
        Private method for running the MCMC simulation
        without interruption.

        number_of_trial_steps : int
           number of trial steps to run without stopping.
        """
        for _ in range(number_of_trial_steps):
            self.do_trial_step()

    def _observe_configuration(self, step):
        """
        Submit current configuration to observers and append observations to
        data container.

        Parameters
        ----------
        step : int
            the current trial step
        """
        row_dict = {}
        new_observations = False
        for observer in self.observers.values():
            if step % observer.interval == 0:
                new_observations = True
                if observer.return_type is dict:
                    for key, value in observer.get_observable(
                            self.calculator.atoms).items():
                        row_dict[key] = value
                else:
                    row_dict[observer.tag] = observer.get_observable(
                        self.calculator.atoms)
        if new_observations:
            self._data_container.append(mctrial=step, record=row_dict)

    @abstractmethod
    def do_trial_step(self):
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
        return random.random()

    @property
    def minimum_observation_interval(self):
        """
        int : minimum number of steps to run
        Monte Carlo simulation without observation interruptions.
        """
        return self._minimum_observation_interval

    def _find_minimum_observation_interval(self):
        """
        Find the greatest common denominator from the observation intervals.
        """

        intervals = [obs.interval for _, obs in self.observers.items()]
        self._minimum_observation_interval = self._get_gcd(intervals)

    def _get_gcd(self, interval_list):
        """
        Find the gcd from the list of intervals.

        interval_list : list of int
        """
        if len(interval_list) == 1:
            return interval_list[0]

        if len(interval_list) > 2:
            gcd_right = gcd(interval_list[-1], interval_list[-2])
            interval_list.pop()
            interval_list.pop()
            interval_list.append(gcd_right)
            return self._get_gcd(interval_list)
        else:
            return gcd(interval_list[0], interval_list[1])

    def attach_observer(self, observer, tag=None):
        """
        Attaches an observer to the ensemble.

        Parameters:
        ----------
        observer : mchammer Observer object
        """
        assert isinstance(observer, BaseObserver), \
            'observer argument invalid; must be child of BaseObserver'

        if tag is not None:
            observer.tag = tag
            self.observers[tag] = observer
        else:
            self.observers[observer.tag] = observer

        if observer.return_type is dict:
            for key in observer:
                self._data_container.add_observable(key, float)
        else:
            self._data_container.add_observable(
                observer.tag, observer.return_type)

        self._find_minimum_observation_interval()

    def reset_data_container(self):
        """Reset the data container and the internal step attribute."""
        self._step = 0
        self._data_container.reset()

    def update_occupations(self, list_of_sites, list_of_elements):
        """
        Update the element occupation of the configuration being sampled.
        This will change the state in both the configuration in the calculator
        and the state of configuration manager.

        parameters
        ----------
        list_of_sites : list of int
            list of indices of the configuration to change.
        list_of_elements : list of int
            list of elements to put on the lattice sites the
            indices refer to.

        raises
        ------
        ValueError : if list_of_sites are not the same length
                     as list_of_elements
        """

        if len(list_of_sites) != len(list_of_elements):
            raise ValueError(
                "List of sites and list of elements are not the same size.")
        self.calculator.update_occupations(list_of_sites, list_of_elements)
        self.configuration.update_occupations(list_of_sites, list_of_elements)
