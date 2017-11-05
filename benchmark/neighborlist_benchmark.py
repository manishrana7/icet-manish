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
from tests.test_manybody_neighborlist import manybodyNeighborlistTester


def benchmark_cpp_nl(structure, cutoff):
    """
    Finds all the neighboring indices up to "order" within the cutoffs using the c++ 
    implemented neighborlist
    """
    
    nl = Neighborlist(cutoff)
    nl.build(structure)


def benchmark_python_nl(atoms, cutoff):
    """
    Finds all the indices up to "order" within the cutoffs using the ase implemented neighborlist
    """
    ase_nl = asenl.NeighborList(len(atoms) * [cutoff / 2.0], skin=1e-8,
                                bothways=True, self_interaction=False)
    ase_nl.update(atoms)


if __name__ == "__main__":
    atoms = bulk("Al").repeat(2)
    structure = structure_from_atoms(atoms)
    cutoff = 10

    t = time.process_time()
    benchmark_cpp_nl(structure, cutoff)
    cpp_elapsed_time = time.process_time() - t
    print("Time taken for cpp neighborlist with  cutoff = {}, len(atoms) = {}".format(
        cutoff, len(atoms)))
    print("{} s".format(cpp_elapsed_time))

    t = time.process_time()
    benchmark_python_nl(atoms,  cutoff)
    py_elapsed_time = time.process_time() - t
    print("Time taken for ASE neigborlist with cutoff = {}, len(atoms) = {}".format(
        cutoff, len(atoms)))
    print("{} s".format(py_elapsed_time))
    print("cpp speedup: {}".format(py_elapsed_time / cpp_elapsed_time))
