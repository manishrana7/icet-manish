from mchammer.calculators.cluster_expansion_calculator import ClusterExpansionCalculator
from icet import ClusterSpace, ClusterExpansion
from ase.build import bulk
import numpy as np
import time


def time_local_energy(calculator, iters=100):
    """ Get time of local energy calculations."""
    indices = [2, 16]
    t0 = time.time()
    calculator.calculate_local_contribution(indices)
    t1 = time.time() - t0
    t1 /= iters
    return t1


def time_total_energy(calculator, iters=100):
    """ Get time of local energy calculations."""
    t0 = time.time()
    calculator.calculate_total()
    t1 = time.time() - t0
    t1 /= iters
    return t1


if __name__ == '__main__':
    atoms = bulk("Al")
    cutoffs = [10, 6, 5]
    elements = ["Al", "Ga"]
    cs = ClusterSpace(atoms, cutoffs, elements)
    parameters = np.array([1.2 for _ in range(len(cs))])
    ce = ClusterExpansion(cs, parameters)
    calculator = ClusterExpansionCalculator(atoms.repeat(10), ce)

    t_local = time_local_energy(calculator)
    t_total = time_total_energy(calculator)
    print("atoms size {}".format(len(calculator.atoms)))
    print("Time taken for local energy {:0.5f}".format(t_local))
    print("Time taken for total energy {:0.5f}".format(t_total))
    print("Speed up for local calc {:0.2f} ".format(t_total/t_local))

    print("Time for calculating one mc step (1000) local energies {}".format(t_local * len(calculator.atoms)))
    print("Time for 1000 mc steps {} minutes".format(t_local * 1000/60* len(calculator.atoms)))
