import timeit
import pyperf
from ase.build import bulk
from icet.core.cluster_space import ClusterSpace as ClusterSpace_cpp

def test():
    atoms = bulk('Al')

    cutoffs = [10, 7, 6]
    chemical_symbols = ['Al', 'Ti']

    cs = ClusterSpace_cpp(atoms, cutoffs, chemical_symbols)  # noqa


if __name__ == '__main__':

#    runner = pyperf.Runner()
#    runner.bench_func('Cluster Space', test)

    print(timeit.timeit("test()",
          setup="from __main__ import test",
          number=100))


