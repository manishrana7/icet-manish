import os
import random
from time import time
from abc import ABC, abstractmethod
from math import gcd

from mchammer.data_container import DataContainer
from mchammer.configuration_manager import ConfigurationManager
from mchammer.observers.base_observer import BaseObserver

import numpy as np

from typing import List, Dict


class BaseEnsemble(ABC):
    """
    Base ensemble abstract class.


    Parameters
    ----------
    calculator : mchammer calculator
        this is the calculator that will be used
        to calculate potential energy differences
        that goes into the Boltzmann factor used
        to calculate the probability ratio between
        two states when doing a trial move.
    atoms : ASE Atoms
        this defines the underlying lattice
        that will be used in the monte carlo
        simulation and also defines the
        initial occupation vector.
    name : str (default BaseEnsemble)
        name of the ensemble which will show up
        in the datacontainer among other things.
    data_container : str
        this defines the filename that the ensemble
        will write the data container to. If the file
        already exists then the data container will be
        read and data will be appended to
        that data container and it will also overwrite
        that data container.
    ensemble_data_write_interval : int
        this sets the interval in which the ensemble
        specific data is observed and saved to the
        data container. This data can be temperature,
        current value of the calculator etc.
    data_container_write_period : float (default np.inf)
        this sets the period in units of seconds
        which the data container should be written to file.
        Writing often to file is both a way to examine
        the currently observed values and also to backup
        your data.
    random_seed : int
        this set the random seed to the random number generator
        that is used in the Monte Carlo trial moves. Setting
        the seed can be used to reproduce earlier
        simulations.

    Attributes
    ----------
    name : str
        name of the ensemble.
    accepted_trials : int
        number of accepted trial steps.
    total_trials : int
        number of total trial steps.
    """

    def __init__(self, calculator=None, atoms=None, name='BaseEnsemble',
                 data_container=None, ensemble_data_write_interval=None,
                 data_container_write_period=np.inf,
                 random_seed=None):

        if calculator is None:
            raise TypeError("Missing required keyword argument: calculator")

        self._calculator = calculator
        self._name = name
        self.data_container_write_period = data_container_write_period
        self.accepted_trials = 0
        self.total_trials = 0
        self._observers = {}
        self._step = 0
        if atoms is None:
            raise TypeError("Missing required keyword argument: atoms")
        if random_seed is None:
            self._random_seed = random.randint(0, 1e16)
        else:
            self._random_seed = random_seed
        random.seed(a=self._random_seed)

        self._data_container_filename = data_container

        if data_container is not None and os.path.isfile(data_container):
            self._data_container = DataContainer.read(data_container)
        else:
            self._data_container = \
                DataContainer(atoms=atoms,
                              ensemble_name=name,
                              random_seed=self._random_seed)

        strict_constraints = self.calculator.occupation_constraints
        sublattices = [[i for i in range(len(self.calculator.atoms))]]
        self.configuration = ConfigurationManager(
            atoms, strict_constraints, sublattices)

        # Handle ensemble data writing
        if ensemble_data_write_interval is None:
            self._ensemble_data_write_interval = len(atoms)
        else:
            self._ensemble_data_write_interval = ensemble_data_write_interval
        self._data_container.add_observable('energy')
        if ensemble_data_write_interval is not np.inf:
            self._find_observer_interval()

    @property
    def structure(self):
        """
        ASE Atoms object :
        The current state of the  structure being sampled in the ensemble.
        """
        return self.configuration.atoms.copy()

    @property
    def step(self) -> int:
        """
        int : current MC trial step
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
    def acceptance_ratio(self) -> float:
        """
        acceptance ratio,
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
        Samples the ensemble for `number_of_trial_steps` steps.

        Parameters
        ----------
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
                first_run_interval = self.observer_interval -\
                    (initial_step -
                     (initial_step // self.observer_interval) *
                        self.observer_interval)
                first_run_interval = min(
                    first_run_interval, number_of_trial_steps)
                self._run(first_run_interval)
                initial_step += first_run_interval
                self._step += first_run_interval

        step = initial_step
        while step < final_step:
            uninterrupted_steps = min(
                self.observer_interval, final_step - step)
            if self._step % self.observer_interval == 0:
                self._observe_configuration(self._step)
            if self._data_container_filename is not None and \
                    time()-last_write_time > self.data_container_write_period:
                self.data_container.write(self._data_container_filename)

            self._run(uninterrupted_steps)
            step += uninterrupted_steps
            self._step += uninterrupted_steps

        # If we end on an observation interval we also observe
        if self._step % self.observer_interval == 0:
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
        Submits current configuration to observers and append observations to
        data container.

        Parameters
        ----------
        step : int
            the current trial step
        """
        row_dict = {}

        # Ensemble specific data
        if step % self._ensemble_data_write_interval == 0:
            ensemble_data = self.get_ensemble_data()
            for key, value in ensemble_data.items():
                row_dict[key] = value

        # Observer data
        for observer in self.observers.values():
            if step % observer.interval == 0:
                if observer.return_type is dict:
                    for key, value in observer.get_observable(
                            self.calculator.atoms).items():
                        row_dict[key] = value
                else:
                    row_dict[observer.tag] = observer.get_observable(
                        self.calculator.atoms)
        if len(row_dict) > 0:
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
        Returns the next random number from the RNG.
        """
        return random.random()

    @property
    def observer_interval(self):
        """
        int : minimum number of steps to run
        Monte Carlo simulation without observation interruptions.
        """
        return self._observer_interval

    def _find_observer_interval(self):
        """
        Find the greatest common denominator from the observation intervals.
        """

        intervals = [obs.interval for obs in self.observers.values()]
        if self._ensemble_data_write_interval is not np.inf:
            intervals.append(self._ensemble_data_write_interval)

        self._observer_interval = self._get_gcd(intervals)

    def _get_gcd(self, interval_list):
        """
        Finds the gcd from the list of intervals.

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

        Parameters
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
            for key in observer.get_keys():
                self._data_container.add_observable(key)
        else:
            self._data_container.add_observable(observer.tag)

        self._find_observer_interval()

    def reset_data_container(self):
        """Resets the data container and the internal step attribute."""
        self._step = 0
        self.total_trials = 0
        self.accepted_trials = 0

        self._data_container.reset()

    def update_occupations(self, list_of_sites, list_of_elements):
        """
        Updates the element occupation of the configuration being sampled.
        This will change the state in both the configuration in the calculator
        and the state of configuration manager.

        Parameters
        ----------
        list_of_sites : list of int
            list of indices of the configuration to change.
        list_of_elements : list of int
            list of elements to put on the lattice sites the
            indices refer to.

        Raises
        ------
        ValueError : if list_of_sites are not the same length
                     as list_of_elements
        """

        if len(list_of_sites) != len(list_of_elements):
            raise ValueError(
                "List of sites and list of elements are not the same size.")
        self.configuration.update_occupations(list_of_sites, list_of_elements)

    def get_property_change(self, indices: List[int],
                            elements: List[int]) -> float:
        """
        Get the property change for a hypothetical change
        of the input indices, elements.

        Parameters
        ----------
        indices : list of int
            refer to indices in the lattice
        elements : list of int
            refer to atomic species

        Returns
        -------
        change in property value : float
        """
        current_elements = [self.configuration.occupations[i] for i in indices]

        current_property = self.calculator.calculate_local_contribution(
            indices, self.configuration.occupations)
        self.update_occupations(list_of_sites=indices,
                                list_of_elements=elements)
        new_property = self.calculator.calculate_local_contribution(
            indices, self.configuration.occupations)
        property_change = new_property - current_property

        # Set elements back to what they were
        self.update_occupations(indices, current_elements)
        return property_change

    def get_ensemble_data(self) -> Dict:
        """
        Get current calculator property.

        Returns
        -------
        dict : ensemble data key pairs

        """
        return {'energy': self.calculator.calculate_total(
            occupations=self.configuration.occupations)}

    def get_random_sublattice_index(self) -> int:
        """
        Get a random sublattice index based
        on the weights of the sublattice.

        Return
        ------
        sublattice_index : int
        """
        return 0
