import numpy as np
import copy
import itertools
from icet.tools.geometry import get_permutation


class Orbit(object):
    """
    Class Orbit

    contains a list of equivalent list of LatticeSites
    contains a sorted Cluster for representation
    Can be compared to other orbits.

    TODO
    ----
    * think about adding __hash__ ?
    * think about overloading orbit + orbit
    Blocked TODO's by Cluster class:
    * geometrical_size() (returns -1 now so we can test it)
    * representative_cluster
    """

    def __init__(self):
        self._equivalent_sites = []
        self._representative_cluster = None
        self.geometrical_size_tolerance = 1e-5
        self._allowed_permutations = []
        self._permutations_to_representative = []

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
    def permutations_to_representative(self):
        """
        Get the list of permutations.
        Where permutations_to_representative[i]
        takes self.equivalent_sites[i] to
        the same order as self.representative_sites.

        Explanation
        -------
        This can be used if you for example want to
        count elements and are interested in difference
        between ABB, BAB, BBA and so on. If you count the
        lattice sites that are permutated according to
        these permutations then you will get the correct
        counts.

        NOTE
        ----
        * This is called _equivalentSitesPermutations
          in C++

        """
        return self._permutations_to_representative

    @permutations_to_representative.setter
    def permutations_to_representative(self, permutations):
        self._permutations_to_representative = permutations

    @property
    def allowed_permutations(self):
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

        NOTE
        ----
        * This is called _allowedSitesPermutations
          in C++
        """
        return self._allowed_permutations

    @allowed_permutations.setter
    def allowed_permutations(self, permutations):
        self._allowed_permutations = permutations

    @property
    def permutated_sites(self):
        """
        Get the equivalent sites but permutated
        to representative site.
        """
        return [self.get_permutated_sites(index) for index in range(len(self))]

    def get_permutated_sites(self, index):
        """
        Return the permutated to representative
        sites of equivalent_sites[index].
        """
        return get_permutation(self.equivalent_sites[index],
                               self.permutations_to_representative[index])

    def get_mc_vectors(self, allowed_components):
        """
        Return the mc vectors for this orbit given the allowed components.
        The mc vectors are returned as a list of tuples

        parameters
       ----------
        allowed_components : list of int
           The allowed components for the lattice sites,
           allowed_components[i] correspond to the lattice site
           self.representative_sites[i].


        """

        assert len(allowed_components) == self.order

        all_possible_mc_vectors = \
            self.get_all_possible_mc_vector_permutations(allowed_components)

        all_possible_mc_vectors.sort()
        mc_vectors = []
        for mc_vector in all_possible_mc_vectors:
            permutated_mc_vector = []
            for allowed_permutation in self.allowed_permutations:
                permutated_mc_vector.append(
                    tuple(get_permutation(mc_vector, allowed_permutation))
                )
            # If this mc vector or any of its allowed permutations
            # exist in mc_vectors append this to the mc vectors
            if set(permutated_mc_vector).isdisjoint(mc_vectors):
                mc_vectors.append(mc_vector)
        return mc_vectors

    def get_all_possible_mc_vector_permutations(self, allowed_components):
        """
       Similar to get all permutations but
       needs to be filtered through the
       number of allowed elements.

       parameters
       ----------
       allowed_components : list of int
           The allowed components for the lattice sites,
           allowed_components[i] correspond to the lattice site
           self.representative_sites[i].

        returns all_mc_vectors : list of tuples of int
        """

        cartesian_lists = []
        for ac in allowed_components:
            cartesian_lists.append(
                [i for i in range(ac - 1)]
            )

        all_mc_vectors = []
        for element in itertools.product(*cartesian_lists):
            all_mc_vectors.append(tuple(element))
        return all_mc_vectors
