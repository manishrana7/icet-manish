from icetdev import *
from icetdev.structure import *
from icetdev.manybodyNeighborlist import *
import numpy.random as random
import numpy as np
from ase import Atoms
#from ase.neighborlist import NeighborList
from ase.build import bulk
import time
import ase.neighborlist as asenl
from tests import manybodyNeighborlistTester


def benchmark_cpp_mbnl(structure, order, cutoff):
    """
    Finds all the indices up to "order" within the cutoffs using the c++ implemented mbnl
    """
    nl = Neighborlist(cutoff)
    nl.build(structure)
    mbnl = ManybodyNeighborlist()
    for i in range(structure.size()):
        mbnl.build(nl, i, order, False)


def benchmark_python_mbnl(atoms, order, cutoff):
    """
    Finds all the indices up to "order" within the cutoffs using the python implemented mbnl
    """
    mbnl_T = manybodyNeighborlistTester.manybodyNeighborlistTester()
    ase_nl = asenl.NeighborList(len(atoms) * [cutoff / 2.0], skin=1e-8,
                                bothways=True, self_interaction=False)
    ase_nl.update(atoms)

    bothways = False
    for i in range(len(atoms)):
        mbnl_T.build(ase_nl, i, order, bothways=bothways)


if __name__ == "__main__":
    atoms = bulk("Al").repeat(2)
    structure = structure_from_atoms(atoms)
    order = 3
    cutoff = 10

    t = time.process_time()
    benchmark_cpp_mbnl(structure, order, cutoff)
    cpp_elapsed_time = time.process_time() - t
    print("Time taken for cpp mbnl with order = {}, cutoff = {}, len(atoms) = {}".format(
        order, cutoff, len(atoms)))
    print("{} s".format(cpp_elapsed_time))

    t = time.process_time()
    benchmark_python_mbnl(atoms, order, cutoff)
    py_elapsed_time = time.process_time() - t
    print("Time taken for python mbnl with order = {}, cutoff = {}, len(atoms) = {}".format(
        order, cutoff, len(atoms)))
    print("{} s".format(py_elapsed_time))
    print("cpp speedup: {}".format(py_elapsed_time / cpp_elapsed_time))
