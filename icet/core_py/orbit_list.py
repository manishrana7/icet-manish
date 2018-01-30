from icet.core_py.permutation_matrix import PermutationMatrix


class OrbitList(object):
    """
    The orbit list object has an internal list of orbits.

    An orbit has a list of equivalent sites with the restriction
    that at least one site is in the primitive cell of the
    primitive structure.

    parameters
    ----------
    atoms : ASE Atoms object
            This atoms object will be used
            to construct a primitive structure
            on which all the lattice sites in the orbits
            are based on.
    cutoffs : list of float
              cutoffs[i] is the cutoff for
              orbits with order i+2.
    verbosity : int
                Set the verbosity for OrbitList and
                all the methods it calls.
    TODO
    ----
    * Write the constructor that should create
      a ManyBodyNeighborList (MBNL) object and a permutation matrix (PM).
      From these two the clusters in the MBNL should find rows in the PM.
      The group of sites made by extracting the columns of the rows make
      up one orbit.
    * The orbit list should be able to sort the orbits to make sure each
      orbit index is the same each time it is constructed with the same
      lattice and cutoffs.
    * Look at the C++ class and mimic useful functionality from there.

    """

    def __init__(self, atoms, cutoffs, verbosity=False):        
        self.permutation_matrix = PermutationMatrix(atoms,max(cutoffs))
        self._primitive_structure = self.permutation_matrix.get_primitive_structure()
        self._orbits = []
        

    def sort(self):
        """
        Sort the list of orbits.
        """
        pass

    @property
    def primitive_structure(self):
        """
        Returns the primitive structure to which the
        lattice sites in the orbits are referenced to.
        """
        return self._primitive_structure().copy()
