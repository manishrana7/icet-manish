from icet.core_py.permutation_matrix import PermutationMatrix

import time
from ase.build import bulk
from icet.core.permutation_map import permutation_matrix_from_atoms
from icet.core.permutation_map import PermutationMap as PermutationMatrix_cpp

from icet.core.orbit_list import __get_lattice_site_permutation_matrix\
    as get_lattice_site_permutation_matrix


def init_permutation_matrix_python(atoms, cutoff):
    """
    Initialize a permutation matrix.
    atoms : ASE Atoms object
    cutofff : float
        defines the neighbor cutoff that will be used
    """

    t = time.process_time()
    pm = PermutationMatrix(atoms, cutoff)
    elapsed_time = time.process_time() - t
    return elapsed_time


def init_permutation_matrix_cpp(atoms, cutoff):
    """
    Initialize a permutation matrix.
    atoms : ASE Atoms object
    cutofff : float
        defines the neighbor cutoff that will be used
    """

    t = time.process_time()
    pm_cpp, prim_structure_cpp, _ = \
        permutation_matrix_from_atoms(
            atoms, cutoff)
    pm_lattice_sites_cpp = \
        get_lattice_site_permutation_matrix(prim_structure_cpp,
                                            pm_cpp,
                                            prune=True)
    elapsed_time = time.process_time() - t
    return elapsed_time


if __name__ == "__main__":

    atoms = bulk("Al")

    cutoff = 15
    python_time = init_permutation_matrix_python(atoms, cutoff)

    cpp_time = init_permutation_matrix_cpp(atoms, cutoff)
    print("Time to initialize pm in pure python with cutoff: {}, {}s ".format(
        cutoff, python_time))
    print("Time to initialize pm in C++/python with cutoff: {}, {}s ".format(
        cutoff, cpp_time))
    print("Speed up in C++ {}".format(python_time/cpp_time))