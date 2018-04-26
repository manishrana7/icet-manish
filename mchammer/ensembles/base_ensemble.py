import random
from time import time
from abc import ABC, abstractmethod
from math import gcd

from mchammer.data_container import DataContainer
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

    @property
    def structure(self):
        """
        ASE Atoms object :
        The current state of the  structure being sampled in the ensemble.
        """
        return self.calculator.atoms.copy()

    @property
    def data_container(self):
        """
        mchammer DataContainer.
        """
        return self._data_container

    @property
    def observers(self):
        """
        dict of mchammer observers.
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

    def run(self, number_of_trial_steps):
        """
        Sample the ensemble for `number_of_trial_steps` steps.

        Parameters:
        -----------

        number_of_trial_steps: int
            number of steps to run in total
        """

        last_backup_time = time()
        for step in range(0, number_of_trial_steps,
                          self.minimum_observation_interval):
            self._run(self.minimum_observation_interval)
            if step % self.minimum_observation_interval == 0:
                self._observe_configuration(step)
            if time() - last_backup_time > self.data_container_write_period \
                    and self._data_container_filename is not None:
                self.data_container.write(self._data_container_filename)

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

    def _observe_configuration(self, step):
        """
        Observe the configuration with the observers that has an interval
        that coincicide with the step.

        parameters
        ----------
        step : int
            the current step
        """
        row_dict = {}
        new_observations = False
        for observer in self.observers.values():
            if step % observer.interval == 0:
                new_observations = True
                if observer.return_type() is dict:
                    for key, value in observer.get_observable(
                            self.calculator.atoms).items():
                        row_dict[key] = value
                else:
                    row_dict[observer.tag] = observer.get_observable(
                        self.calculator.atoms)
        if new_observations:
            self._data_container.append(mctrial=step, record=row_dict)

    def _find_minimum_observation_interval(self):
        """
        Find the greatest common denominator from the observation intervals.
        """

        intervals = [obs.interval for _, obs in self.observers.items()]
        if len(intervals) == 1:
            self._minimum_observation_interval = intervals[0]
        else:
            self._minimum_observation_interval = self._get_gcd(intervals)

    def _get_gcd(self, interval_list):
        """
        Find the gcd from the list of intervals.

        interval_list : list of int
        """

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

        if observer.return_type() is dict:
            for key in observer:
                self._data_container.add_observable(key, float)
        else:
            self._data_container.add_observable(
                observer.tag, observer.return_type())

        self._find_minimum_observation_interval()
