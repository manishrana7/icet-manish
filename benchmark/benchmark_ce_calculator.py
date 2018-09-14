from mchammer.calculators.cluster_expansion_calculator import ClusterExpansionCalculator
from icet import ClusterSpace, ClusterExpansion
from ase.build import bulk
import numpy as np
import time


def time_local_energy(calculator, iters):
    """ Get time of local energy calculations."""
    indices = [2]
    occupations = calculator.atoms.numbers
    t0 = time.time()
    for _ in range(iters):
        calculator.calculate_local_contribution(local_indices=indices, occupations=occupations)
    t1 = time.time() - t0
    t1 /= iters
    return t1


def time_total_energy(calculator, iters):
    """ Get time of local energy calculations."""
    occupations = calculator.atoms.numbers
    t0 = time.time()
    for _ in range(iters):
        calculator.calculate_total(occupations=occupations)
    t1 = time.time() - t0
    t1 /= iters
    return t1


if __name__ == '__main__':
    iters=100
    atoms = bulk("Al")
    cutoffs = [10, 6, 5]
    elements = ["Al", "Ga"]
    cs = ClusterSpace(atoms, cutoffs, elements)
    atoms = atoms.repeat(10)
    parameters = np.array([1.2 for _ in range(len(cs))])
    ce = ClusterExpansion(cs, parameters)
    print("Beging construct ce calc")
    t0 = time.time()
    calculator = ClusterExpansionCalculator(atoms, ce)
    print("done construct ce calc. Time: {:0.5f}s".format(time.time()-t0))
    t_local = time_local_energy(calculator, iters)
    t_total = time_total_energy(calculator, iters)
    print("atoms size {}".format(len(calculator.atoms)))
    print("Time taken for local energy {:0.5f}".format(t_local))
    print("Time taken for total energy {:0.5f}".format(t_total))
    print("Speed up for local calc {:0.2f} ".format(t_total/t_local))

    print("Time for calculating {} mc step ({}) local energies {}s".format(iters/len(atoms),len(atoms),t_local * len(calculator.atoms)))
    print("Time for {} mc steps {} minutes".format(iters/len(atoms), t_local /60 * iters/len(calculator.atoms)))
