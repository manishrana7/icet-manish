import numpy as np
import copy


class Orbit(object):
    """
    Class Orbit

    contains a list of equivalent list of LatticeSites
    contains a sorted Cluster for representation
    Can be compared to other orbits.

    TODO
    ----
    * write constructor
    * Add functions
        * getPermutatedEquivalentSites
        * getSitesWithPermutation
        * get_mc_vectors
    * properties

    * think about adding __hash__ ?
    * think about overloading orbit + orbit
    Blocked TODO's by Cluster class:
    * geometrical_size()
    * representative_cluster
    """

    def __init__(self):
        self._equivalent_sites = []
        self._representative_cluster = None
        self.geometrical_size_tolerance = 1e-5
        self._equivalent_permutations = []

    @property
    def equivalent_sites(self):
        """
        List of equivalent Lattice Sites
        """
        return self._equivalent_sites

    @equivalent_sites.setter
    def equivalent_sites(self, sites):
        self._equivalent_sites = sites

    @property
    def representative_cluster(self):
        """
        The representative cluster
        represents the geometrical
        version of what this orbit is.
        """
        return self._representative_cluster

    @property
    def representative_sites(self):
        """
        The representative sites
        is a list of lattice sites
        that are uniquely picked out
        for this orbit which can be
        used to represent and distinguish
        between orbits.
        """
        if len(self.equivalent_sites) > 0:
            return self.equivalent_sites[0]
        else:
            raise IndexError('Equivalent sites is empty')

    @property
    def order(self):
        """
        Returns the order of the orbit.
        The order is the same as the number
        of bodies in the representative cluster
        or the number of lattice sites per element
        in equivalent_sites.
        """
        return len(self.representative_sites)

    def __len__(self):
        """
        Lenght of orbit is defined as the
        number of equivalent sites.
        """
        return len(self.equivalent_sites)

    @property
    def geometrical_size(self):
        """
        Returns the geometrical size of the
        representative cluster.
        TODO
        ----
        * Implement this when cluster is available.
        """
        return -1

    def sort(self):
        """
        Sort the equivalent sites list.
        """
        self._equivalent_sites.sort()

    def __eq__(self, other):
        """
        Test equivalence between two orbits.
        The orbits will be marked as the same
        if all the equivalent_sites are identical.
        """
        return self.equivalent_sites == other.equivalent_sites

    def __lt__(self, other):
        """
        less than operator for sorting in containers.
        It will compare properties in this order:
        * order
        * geometrical size
        * len(equivalent_sites)
        * finally the equivalent sites in lexicographical order
        """
        if self.order != other.order:
            return self.order < other.order

        if np.abs(self.geometrical_size -
                  other.geometrical_size) > self.geometrical_size_tolerance:
            return self.geometrical_size < other.geometrical_size

        if len(self) != len(other):
            return len(self) < len(other.orbit)

        return self.equivalent_sites < other.equivalent_sites

    def __add__(self, other):
        """
        Define the add operator.
        Allowed values:
        * type ndarray with shape(3,)
        """
        if not isinstance(other, type(np.array([1, 1, 1]))) or len(other) != 3:
            raise TypeError("Adding orbit with {}".format(other))

        orbit = Orbit()
        orbit.equivalent_sites = copy.deepcopy(self.equivalent_sites)
        for sites in orbit.equivalent_sites:
            for site in sites:
                site.unitcell_offset += other
        return orbit

    @property
    def equivalent_permutations(self):
        """
        Get the list of equivalent permutations
        for this orbit.

        Explanation
        -------
        If this orbit is a triplet
        and the permutation [0,2,1]
        exists this means that
        The lattice sites [s1, s2, s3]
        are equivalent to [s1, s3, s2]
        This will have the effect that
        for a ternary CE the cluster
        functions (0,1,0) will not
        be considered since it is
        equivalent to (0,0,1).
        """
        return self._equivalent_permutations

    @equivalent_permutations.setter
    def equivalent_permutations(self, permutations):
        self._equivalent_permutations = permutations
