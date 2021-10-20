import numpy as np
import time
from ase.build import bulk

from icet import ClusterSpace, ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator


def time_local_energy(calculator, occupations, iters):
    """ Returns timing for energy change calculations. """
    t0 = time.time()
    for _ in range(iters):
        calculator._calculate_partial_change(occupations=occupations,
                                             flip_index=2,
                                             new_occupation=31)
    t1 = time.time() - t0
    return t1


def time_total_energy(calculator, occupations, iters):
    """ Returns timing for total energy calculations. """
    t0 = time.time()
    for _ in range(iters):
        calculator.calculate_total(occupations=occupations)
    t1 = time.time() - t0
    return t1


def print_timing_ratios(structure, local_iters, total_iters, sizes, cutoffs):
    """ Prints timing ratios between local and total energy calculations. """
    print('# Col 1: Structure side length')
    print('# Col 2: Number of atoms in structure')
    print('# Col 3: Number of iterations for change calculation')
    print('# Col 4: Number of iterations for total calculation')
    print('# Col 5: Initialization time (s)')
    print('# Col 6: Total time (all iterations) for change calculation (s)')
    print('# Col 7: Total time (all iterations) for total calculation (s)')
    print('# Col 8: Time for one total calculation divided by one change calculation')

    cs = ClusterSpace(structure, cutoffs, chemical_symbols=['Al', 'Ga'])
    parameters = np.array([1.2 for _ in range(len(cs))])
    ce = ClusterExpansion(cs, parameters)
    for size in sizes:
        structure_cpy = structure.repeat(size)
        occupations = structure_cpy.get_atomic_numbers()
        t0 = time.time()
        calculator = ClusterExpansionCalculator(structure_cpy, ce)
        time_ce_init = time.time() - t0
        t_local = time_local_energy(calculator, occupations, local_iters)
        t_total = time_total_energy(calculator, occupations, total_iters)
        print(f'{size:3d} {len(structure_cpy):5d} {local_iters:5d} {total_iters:2d} '
              f'{time_ce_init:10.6f} {t_local:10.6f} {t_total:10.6f} '
              f'{(t_total / total_iters) / (t_local / local_iters):12.4f}')


if __name__ == '__main__':

    local_iters = 10000
    total_iters = 100
    structure = bulk('Al')
    cutoffs = [10, 6]
    chemical_symbols = ['Al', 'Ga']
    sizes = [4, 6, 8, 10, 14, 16, 20]
    print_timing_ratios(structure, local_iters, total_iters, sizes, cutoffs)
