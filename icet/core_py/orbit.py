

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
        * Add equivalent sites
        * __len__
        * geometrical_size()        
        * getPermutatedEquivalentSites
        * getSitesWithPermutation        
        * sort
        * __eq__
        * __lt__
        * __hash__ ?
        * orbit + offset
        * orbit + orbit
        * get_mc_vectors
    * properties       
       * order       
    """

    def __init__(self):
        self._equivalent_sites = []
        self._representative_cluster = None

    @property
    def equivalent_sites(self):
        """
        List of equivalent Lattice Site's        
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
        if len(self.equivalent_sites)>0:
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