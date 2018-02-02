import numpy as np

from icet.core_py.permutation_matrix import PermutationMatrix
from icet.core_py.orbit import Orbit
from icet.core_py.many_body_neighbor_list import ManyBodyNeighborList
from icet.core_py.lattice_site import LatticeSite



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
        self.column1 = self.permutation_matrix.column1

        self._primitive_structure = self.permutation_matrix.primitive_structure
        mbnl = ManyBodyNeighborList(self.primitive_structure, cutoffs)
        self.taken_rows = set()
        self._orbits = []
        for index in range(len(self.primitive_structure)):
            mb_neigbhors_index = mbnl.build(index)
            for compressed_sites in mb_neigbhors_index:
                for sites in mbnl.unzip(compressed_sites):
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
        return self._primitive_structure.copy()

    @property
    def permutation_matrix(self):
        """
        Return icet PermutationMatrix object.
        """
        return self._permutation_matrix

    @property
    def orbits(self):
        """
        Return the internal list of orbits
        """
        return self._orbits

    def is_new_orbit(self, sites):
        """
        Checks if this sites has been added into
        an orbit allready.
        """
        if len(sites) == 0:
            raise RuntimeError("sites is empty in is new orbit")
        translated_eq_sites = self.get_all_translated_sites(sites)
        if len(translated_eq_sites) == 0:
            raise RuntimeError("translated_eq_sites is empty in is new orbit")

        sites_indices_match = self.get_matches_in_pm(translated_eq_sites)
        if not self.is_rows_taken(sites_indices_match[0][1]):
            return True
        return False

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
        # if len(sites) ==1:
        #     return orbit
        assert len(rows) == len(sites)
        assert isinstance(rows, list)
        assert isinstance(rows[0], list)
        assert isinstance(rows[0][0], LatticeSite), "{} != LatticeSite".format(
            type(rows[0][0]))

        for i in range(len(sites)):
            assert len(rows[i]) == len(rows[0])

        for i in range(len(rows[0])):
            eq_sites = [row[i] for row in rows]
            assert len(eq_sites) == len(sites), "{} != {}".format(
                len(eq_sites), len(sites))
            
            translated_eq_sites = self.get_all_translated_sites(eq_sites)
            sites_indices_match = self.get_matches_in_pm(translated_eq_sites)
            if not self.is_rows_taken(sites_indices_match[0][1]):
                orbit.equivalent_sites.append(eq_sites)                
                for site_index in sites_indices_match:
                    self.take_row(site_index[1])

        return orbit

    def get_rows(self, sites):
        """
        Returns rows of the permutation matrix
        corresponding to the sites in the sites.

        Parameters
        ----------

        sites : list of icet Lattice Site objects        

        """
        row_indices = self.get_row_indices(sites)

        return [self.permutation_matrix.pm_lattice_sites[index] for index in row_indices]

    def get_row_indices(self, sites):
        """
        Return the indices for the rows that the sites
        in the sites matches to.

        TODO : think if this should be sorted
        TODO : Should this be a tuple for easier hashing? Yes!
        """
        indices = [None] * len(sites)

        for i, col_site in enumerate(self.column1):
            for j, site in enumerate(sites):
                if site == col_site:
                    indices[j] = i

        for index in indices:
            if index == None:
                raise RuntimeError("index not found for sites")
        return tuple(sorted(indices))  # TODO check if this should be done elsewhere
        # return [self.column1.index(site) for site in sites]

    def get_all_translated_sites(self, sites):
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
            if not np.allclose(site.unitcell_offset, np.array([0., 0., 0.])):
                translated_sites.append(
                    [ls - site.unitcell_offset for ls in sites])

        if len(translated_sites) == 0:
            return [sites]
        return translated_sites

    @property
    def permutation_matrix(self):
        """
        Returns the core_py permutation matrix object.
        """
        return self._permutation_matrix

    def __len__(self):
        """
        Lenght of an orbit list is number of orbits in orbit list.
        """
        return len(self._orbits)

    def get_matches_in_pm(self, list_of_sites):
        """
        Returns a list of tuple of lattice sites and
        matched rows in pm of those that has a
        match in column 1.

        Parameters
        ----------
        list_of_sites : list of list of LatticeSite

        Return
        ------
        matched_sites : the elements in sites that had
        a match in column 1
        """
        if len(list_of_sites) == 0:
            raise RuntimeError("List of sites empty")
        matched_sites = []
        for sites in list_of_sites:
            try:
                rows = self.get_row_indices(sites)
                matched_sites.append(tuple((sites, rows)))
            except RuntimeError as e:
                if "index not found for sites" in e:
                    continue
                else:
                    raise RuntimeError(e)
        if len(matched_sites) > 0:
            return matched_sites
        else:
            raise RuntimeError("Did not find any of the "
                               "translated sites in col1 "
                               "of permutation matrix in "
                               "function get_matches_in_pm"
                               "in orbit list")

    def __str__(self):
        nice_str = ''
        for i, orbit in enumerate(self.orbits):
            nice_str += "orbit {} - Multiplicity {} '\n'".format(i,len(orbit))
        return nice_str

    def is_rows_taken(self, rows):
        """
        Checks if these particual rows in the
        permutation matrix has been used for
        constructing an orbit.

        parameters
        ---------
        rows : list of ints
            Refers to row indices of
            the permutation matrix
        """        
        if rows in self.taken_rows:
            return True
        return False

    def take_row(self, row):
        """
        Add this row to the list of taken rows
        in the permtuation matrix.

        row : list (tuple) of int        
        """
        self.taken_rows.add(row)
