from icet.core.cluster_space import ClusterSpace as ClusterSpace_cpp
from icet.core_py.cluster_space import ClusterSpace as ClusterSpace_py

from ase.build import bulk

import time


def setup_cluster_space_py(atoms, cutoffs, elements):
    """
    Initialize the python version of the cluster space object.

    parameters
    ----------
    atoms : Ase Atoms object
    cutoffs : list of floats
    elements : list of str
    """

    t = time.process_time()
    cs = ClusterSpace_py(atoms, cutoffs, elements)  # noqa
    elapsed_time = time.process_time() - t
    return elapsed_time


def setup_cluster_space_cpp(atoms, cutoffs, elements):
    """
    Initialize the C++ version of the cluster space object.

    parameters
    ----------
    atoms : Ase Atoms object
    cutoffs : list of floats
    elements : list of str
    """

    t = time.process_time()
    cs = ClusterSpace_cpp(atoms, cutoffs, elements)  # noqa
    elapsed_time = time.process_time() - t
    return elapsed_time


if __name__ == "__main__":

    atoms = bulk("Al")

    cutoff = [15, 10]
    elements = ["Al", "Ti"]
    python_time = setup_cluster_space_py(atoms, cutoff, elements)

    cpp_time = setup_cluster_space_cpp(atoms, cutoff, elements)
    print("Time to initialize ClusterSpace in pure "
          "python with cutoff: {}, {}s ".format(
              cutoff, python_time))
    print("Time to initialize ClusterSpace"
          " in C++/python with cutoff: {}, {}s ".format(
              cutoff, cpp_time))
    print("Speed up in C++ {}".format(python_time / cpp_time))
