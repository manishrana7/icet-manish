
from icet.core_py.orbit_list import OrbitList
from ase.build import bulk
import time

def time_orbit_list(atoms, cutoffs):
    """
    Return timing on creating orbitlist    
    """
    t = time.process_time()
    orbit = OrbitList(atoms, cutoffs)
    elapsed_time = time.process_time() - t

    return elapsed_time




if __name__ == "__main__":

    atoms = bulk("Al")

    cutoff = [10, 7, 6]
    python_time = time_orbit_list(atoms, cutoff)

    print("Time to initialize orbitlist in pure "
          "python with cutoff: {}, {}s ".format(
              cutoff, python_time))
    # print("Time to initialize ClusterSpace"
    #       " in C++/python with cutoff: {}, {}s ".format(
    #           cutoff, cpp_time))
    # print("Speed up in C++ {}".format(python_time / cpp_time))
