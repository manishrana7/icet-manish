from icetdev import *
from icetdev.structure import *
from icetdev.manybodyNeighborlist import *
import numpy.random as random
import numpy as np
from ase import Atoms
#from ase.neighborlist import NeighborList
from ase.build import bulk
import time


def benchmark_cpp_mbnl(structure, order, cutoff):
    nl = Neighborlist(cutoff)
    nl.build(structure)
    mbnl = ManybodyNeighborlist()
    for i in range(structure.size()):
        mbnl.build(nl, i, order, False)


if __name__ == "__main__":
    atoms = bulk("Al").repeat(2)
    structure = structure_from_atoms(atoms)
    order = 5
    cutoff = 10
    t = time.process_time()
    benchmark_cpp_mbnl(structure, order, cutoff)
    elapsed_time = time.process_time() - t
    print("Time taken for cpp mbnl with order = {}, cutoff = {}, len(atoms) = {}".format(
        order, cutoff, len(atoms)))
    print("{} s".format(elapsed_time))
