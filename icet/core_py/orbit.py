

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
       * representative cluster
       * Representative sites
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
        
