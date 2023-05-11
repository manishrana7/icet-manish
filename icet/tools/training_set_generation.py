import numpy as np
from typing import List, Tuple
from mchammer.ensembles.canonical_annealing import available_cooling_functions
from icet import (ClusterSpace, StructureContainer)
from ase import Atoms
from icet.input_output.logging_tools import logger
logger = logger.getChild('training_set_generation')


def _covariance(X: np.ndarray) -> np.ndarray:
    """
    Get the norm of the covariance
    between parameters.

    Ignore diagonal since self correlation is
    not important

    Parameters
    ----------
    X
        current fit matrix
    """
    cov = X.T @ X
    cov = np.tril(cov, k=-1)
    return np.linalg.norm(cov)


def _get_fit_matrix(structure_container: StructureContainer,
                    new_inds: np.ndarray,
                    n_base_structures: int) -> np.ndarray:
    """
    Get the current fit matrix

    Parameters
    ----------
    structure_container
        The structure container
    new_inds
        The part of the structure container that contains the
        new structures to be added
    n_base_structures
        Number of structures in the base
    """
    base_inds = np.array([f for f in range(n_base_structures)], dtype=int)
    inds = np.append(base_inds, n_base_structures + new_inds)
    (fit_matrix, _) = structure_container.get_fit_data(inds)
    return fit_matrix


def _do_swap(inds: np.ndarray,
             n_structures_to_add: int,
             n_mcmc_structures: int) -> np.ndarray:
    """
    Update indices to be used as training data

    Parameters
    ----------
    inds
        The current indicies that are used
    n_structures_to_add
        Total number of structures to add to the current base structures
    n_mcmc_structures
        The size of the pool of potential candidate structures
    """
    # Get index we swap out
    _inds = inds.copy()
    swap_out = np.random.choice(range(inds.size))

    # Get index from the pool that are not currently in inds
    inds_pool = np.array([range(n_mcmc_structures)])
    inds_pool = np.setdiff1d(inds_pool, inds, assume_unique=True)

    # Get index of structure to swap in
    swap_in = np.random.choice(inds_pool)

    # Do the swap
    _inds[swap_out] = swap_in
    return _inds


def structure_annealing(
        cluster_space: ClusterSpace,
        monte_carlo_structures: List[Atoms],
        n_structures_to_add: int,
        n_steps: int,
        base_structures: List[Atoms] = None,
        cooling_start: float = 5,
        cooling_stop: float = 0.001,
        metric_function: str = 'condition number',
        cooling_function: str = 'exponential',
        start_inds: List[int] = None)\
            -> Tuple[List[int], List[float]]:
    """
    Given a cluster space, a base pool of structures and new pool of structures
    use an Monte Carlo inspired annealing method to find a good structure pool
    for training.

    The function returns the indices of the optimal structures in the
    monte_carlo_structures pool and a list of accepted metric values

    Parameters
    ----------
    cluster_space
        A cluster space defining the lattice to be occupied
    monte_carlo_structures
        A list of candidate training structures
    n_structures_to_add
        How many of the structures in the monte_carlo_structures pool that
        should be kept for training
    n_steps
        Number of steps in the annealing algorithm
    base_structures
        A list of structures that is already in your training pool
        can be None if you don't have any structures yet.
    cooling_start
        Initial value of the cooling_function
    cooling_stop
        Last value of the cooling_function
    cooling_function : str/function
        Artificial number that rescales the difference between the metric value
        between two iterations.
        Available options are `linear` or `exponential`
    metric_function : str/function
        Metric to measure how good your current structure pool is.
        Available options are `condition number` or `covariance`
    start_inds
        Picks out the starting structure from the monte_carlo_structures pool.
        Can be used if you want to continue from an old run for example.

    Example
    ----------
        >>> prim = bulk('Au', a=4.0)
        >>> cs = ClusterSpace(prim, [6.0], [['Au', 'Pd']])
        >>> target_concentrations = []
        >>> for i in range(1, 10):
        >>>     target_concentrations.append({'Au': i / 10, 'Pd': (10 - i) / 10})
        >>> structure_pool = []
        >>> for _ in range(int(5e2)):
        >>>     structure = prim.copy().repeat(10)
        >>>     indx = np.random.choice(range(9))
        >>>     concentration = target_concentrations[indx]
        >>>     occupy_structure_randomly(structure, cs,
        >>>                               concentration)
        >>>     structure_pool.append(structure)
        >>> start_inds = [f for f in range(10)]
        >>> (inds, traj) = structure_annealing(cs,
        >>>                                    structure_pool,
        >>>                                    10,
        >>>                                    int(1e3),
        >>>                                    start_inds=start_inds)
        >>> training_structures = []
        >>> for ind in inds:
        >>>     training_structures.append(structure_pool[ind]))

    """
    if base_structures is None:
        base_structures = []

    # setup metric function
    if isinstance(metric_function, str):
        available = sorted(available_metric_functions.keys())
        if metric_function not in available:
            raise ValueError(
                'Select from the available metric_functions {}'.format(available))
        _metric_function = available_metric_functions[metric_function]
    elif callable(metric_function):
        _metric_function = metric_function
    else:
        raise TypeError('metric_function must be either str or a function')

    # setup cooling function
    if isinstance(cooling_function, str):
        available = sorted(available_cooling_functions.keys())
        if cooling_function not in available:
            raise ValueError(
                'Select from the available cooling_functions {}'.format(available))
        _cooling_function = available_cooling_functions[cooling_function]
    elif callable(cooling_function):
        _cooling_function = cooling_function
    else:
        raise TypeError('cooling_function must be either str or a function')

    # setup cluster vectors
    structure_container = StructureContainer(cluster_space)
    for structure in base_structures:
        structure_container.add_structure(structure, properties={'energy': 0})
    for structure in monte_carlo_structures:
        structure_container.add_structure(structure, properties={'energy': 0})

    # get number of structures in monte_carlo_structures
    n_mcmc_structures = len(monte_carlo_structures)
    n_base_structures = len(base_structures)

    # randomly chose starting structure unless user want specific indices
    if not start_inds:
        inds = np.random.choice(range(len(monte_carlo_structures)), size=n_structures_to_add,
                                replace=False)
    else:
        inds = np.array(start_inds)

    # get initial fitting_matrix, A in Ax = y
    fit_matrix = _get_fit_matrix(structure_container, inds, n_base_structures)

    # get metric of fitting_matrix
    metric_val = _metric_function(fit_matrix)
    metric_traj = [metric_val]
    for n in range(n_steps):
        # get current artificial cooling
        T = _cooling_function(n, cooling_start, cooling_stop, n_steps)

        # do swaps from pool
        new_inds = _do_swap(inds, n_structures_to_add, n_mcmc_structures)
        new_fit_matrix = _get_fit_matrix(structure_container, new_inds, n_base_structures)

        # get new metric
        metric_new = _metric_function(new_fit_matrix)
        if -(metric_new - metric_val) / T > np.log(np.random.uniform()):
            # if accepted update data
            metric_val = metric_new
            inds = new_inds
            metric_traj.append(metric_val)

        if (n % 100):
            logger.info(f'step {n:6d}, T {T:8.5f} , current metric {metric_val:8.5f}')

    return (inds, metric_traj)


available_metric_functions = {'covariance': _covariance,
                              'condition number': np.linalg.cond}
