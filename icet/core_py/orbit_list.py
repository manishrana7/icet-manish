from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core_py.orbit import Orbit
from icet.core_py.many_body_neighbor_list import ManyBodyNeighborList


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
      From these two the sites in the MBNL should find rows in the PM.
      The group of sites made by extracting the columns of the rows make
      up one orbit.
    * The orbit list should be able to sort the orbits to make sure each
      orbit index is the same each time it is constructed with the same
      lattice and cutoffs.
    * Look at the C++ class and mimic useful functionality from there.

    * Add cluster to a newly created orbit
    """

    def __init__(self, atoms, cutoffs, verbosity=False):
        self._permutation_matrix = PermutationMatrix(atoms, max(cutoffs))
        self._column1 = self.permutation_matrix.column1
        self._primitive_structure = self.permutation_matrix.get_primitive_structure()
        mbnl = ManyBodyNeighborList(self.primitive_structure, cutoffs)

        self._orbits = []
        for index in range(len(self.primitive_structure)):
            mb_neigbhors_index = mbnl.build(index)
            for sites in mbnl.unzip(mb_neigbhors_index):
                if self.is_new_orbit(sites):
                    orbit = self.make_orbit(sites)
                    self._orbits.append(orbit)

    def sort(self):
        """
        Sort the list of orbits.
        """
        self._orbits.sort()

    @property
    def primitive_structure(self):
        """
        Returns the primitive structure to which the
        lattice sites in the orbits are referenced to.
        """
        return self._primitive_structure().copy()

    @property
    def permutation_matrix(self):
        """
        Return icet PermutationMatrix object.
        """
        return self._permutation_matrix

    def is_new_orbit(self, sites):
        """
        Checks if this sites has been added into
        an orbit allready.
        """
        pass

    def make_orbit(self, sites):
        """
        Creates a new orbit from the sites.
        It will find the rows of the first column in
        the permutation matrix that correspond
        to the sites in the sites and add all
        equivalent sites in an Orbit by traversing the
        columns in the permutation matrix

        Parameters
        ----------
        sites : list of icet Lattice Site objects

        returns orbit : icet Orbit object

        TODO:
        * add cluster (geometrical version) to the orbit
        """
        orbit = Orbit()

        # TODO check sorted rows?
        rows = self.get_rows(sites)

        for eq_sites in zip(rows):
            translated_eq_sites = self.get_all_translated_sites(eq_sites)

    def get_rows(self, sites):
        """
        Returns rows of the permutation matrix
        corresponding to the sites in the sites.

        Parameters
        ----------

        sites : list of icet Lattice Site objects        

        """
        row_indices = self.get_row_indices(sites)

        return [self.permutation_matrix.pm_lattice_site[index] for index in row_indices]

    def get_row_indices(self, sites):
        """
        Return the indices for the rows that the sites
        in the sites matches to.

        TODO : think if this should be sorted
        TODO : Should this be a tuple for easier hashing?
        """
        return [self.column1.index(site) for site in sites]

    def get_all_translated_sites(self, sites)
    """
    Construct a list of lists of sites. 
    Will for each site that has
    a unitcell offset !=[0,0,0] translate all sites so
    *that* site is in [0,0,0] and the others sides just
    go along for the ride.

    The resulted list are all equivalent sites that
    the permutation matrix doesn't really catch on.

    parameters
    ----------
    sites : list of icet Lattice Sites object
    """
    translated_sites = []
    for site in sites:
        if site.unitcell_offset != [0, 0, 0]:
            translated_sites.append(
                [ls + site.unitcell_offset for ls in sites])

    return translated_sites
