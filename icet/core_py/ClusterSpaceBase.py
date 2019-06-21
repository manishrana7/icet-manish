import copy
import numpy as np

from _icet import _OrbitList
from _icet import ClusterCounts

from icet.core.local_orbit_list_generator import LocalOrbitListGenerator
#from icet.core.cluster_counts import ClusterCounts
from icet.core.structure import Structure
from icet.core_py.PeriodicTable import PeriodicTable

class ClusterSpaceBase(object):
    """
    @brief This class handles the cluster space.
    @details It provides functionality for setting up a cluster space, calculating cluster vectors as well as retrieving various types of associated information.
    """
    def __init__(self, chemicalSymbols: list, orbitList: _OrbitList):
        """
        @details This constructor initializes a ClusterSpace object.
        @param chemicalSymbols vector of allowed chemical symbol for each site
        @param orbitList list of orbits for the primitive structure
        """
        self._orbitList = orbitList
        self._chemicalSymbols = chemicalSymbols
        self._clusterSpaceInfo = None
        self._speciesMap = []

        self._numberOfAllowedSpeciesPerSite = \
            [len(self._chemicalSymbols[i]) for i in range(len(self._chemicalSymbols))]

        self._primitiveStructure = self._orbitList.get_primitive_structure()
        self._primitiveStructure.set_number_of_allowed_species(self._numberOfAllowedSpeciesPerSite)

        # Set up a map between chemical species and the internal species enumeration scheme.
        for i in range(len(self._primitiveStructure)):
            speciesMap = {}
            species = []
            for el in chemicalSymbols[0]:
                species.append(PeriodicTable.strInt[el])
            species.sort()

            for i in range(len(species)):
                speciesMap[species[i]] = i

            self._speciesMap.append(speciesMap)

        self.precomputeMultiComponentVectors()


    def __len__(self) -> int:
        return len(self._orbitList)

    # Returns the cutoff for each order.
    def getCutoffs(self):
        return self._clusterCutoffs

    # Returns the primitive structure.
    def getPrimitiveStructure(self):
        return self._primitiveStructure

    # Returns the orbit for an index.
    def getOrbit(self, index:int):
        return self._orbitList.get_orbit(index)

    def getClusterVector(self, structure: Structure):
        """
        @details This function calculates and then returns the cluster vector for the
        input structure in the cluster space.

        The first element in the cluster vector will always be one (1) corresponding to
        the zerolet. The remaining elements of the cluster vector represent averages
        over orbits (symmetry equivalent clusters) of increasing order and size.

        @param structure input configuration
        """

        # Don't sort clusters
        orderIntact = True
        # Count the clusters in the orbit with the same order as the prototype cluster
        permuteSites = True

        localOrbitListGenerator = LocalOrbitListGenerator(self._orbitList, structure)
        uniqueOffsets = localOrbitListGenerator.get_number_of_unique_offsets()
        clusterCounts = ClusterCounts()

        # Create local orbit lists and count associated clusters.
        for i in range(uniqueOffsets):
            localOrbitList = localOrbitListGenerator.generate_local_orbit_list(i)
            clusterCounts.count_orbit_list(structure, localOrbitList, orderIntact, permuteSites)

        # Check that the number of unique offsets equals the number of unit cells in the supercell.
        numberOfUnitcellRepetitions = len(structure)/len(self._primitiveStructure)

        if uniqueOffsets != numberOfUnitcellRepetitions:
            msg = "The number of unique offsets does not match the number of primitive units in the input structure: "
            msg += str(uniqueOffsets) + " != " + str(numberOfUnitcellRepetitions)
            raise RuntimeError(msg)

        # Get the cluster -> cluster counts map
        clusterMap = clusterCounts.get_cluster_counts()

        # Initialize cluster vector and insert zerolet.
        clusterVector = [1]

        # Loop over orbits.
        for i_orbit in range(len(self._orbitList)):
            representativeCluster = self._orbitList.orbits[i_orbit].get_representative_cluster()
            # @todo This is necessary. Side effects need to carefully evaluated.
            # Ideally move this task elsewhere as it is repeated for every structure,
            # for which the cluster vector is computed.
            representativeCluster.tag = i_orbit

            numberOfAllowedSpeciesBySite = []

            try:
                numberOfAllowedSpeciesBySite =\
                    self.getNumberOfAllowedSpeciesBySite(self._primitiveStructure,
                                                    self._orbitList.orbits[i_orbit].get_representative_sites());
            except Exception as e:
                msg = "Failed retrieving the number of allowed species in ClusterSpace.getClusterVector\n"
                msg += str(e)
                raise RuntimeError(msg)

            # Jump to the next orbit if any of the sites in the representative
            # cluster are inactive (i.e. the number of allowed species on this
            # site is less than 2).
            if any([k<2 for k in numberOfAllowedSpeciesBySite]):
                continue

            representativeSites = self._orbitList.orbits[i_orbit].get_representative_sites()
            representativeSitesIndices = [site.index for site in representativeSites]

            # First we obtain the multi-component vectors for this orbit, i.e. a vector
            # of vectors of int (where the int represents a cluster function index).
            # Example 1: For an AB alloy we obtain [0, 0] and [0, 0, 0] for pair and triplet terms, respectively.
            # Example 2: For an ABC alloy we obtain [0, 0], [0, 1], [1, 1] for pairs and similarly for triplets.
            # Depending on the symmetry of the cluster one might also obtain [1, 0]
            # (e.g., in a clathrate or for some clusters on a HCP lattice).
            multiComponentVectors = self._multiComponentVectors[i_orbit]

            # @todo Make getMultiComponentVectorPermutations take an Orbit rather than an index. Then swap the loop over int for a loop over Orbit above.
            sitePermutations = self._sitePermutations[i_orbit]

            currentMultiComponentVectorIndex = 0

            for multiComponentVector in multiComponentVectors:
                clusterVectorElement = 0
                multiplicity = 0

                # speciesCountPair -> {('Ag',): 1}
                speciesCountPair = clusterMap[representativeCluster]
                # convert symbols to atomic numbers
                speciesSymbols = list(speciesCountPair.keys())[0]
                species_an = [PeriodicTable.strInt[s] for s in speciesSymbols]

                speciesValue = list(speciesCountPair.values())[0]

                for perm in sitePermutations[currentMultiComponentVectorIndex]:
                    permutedMultiComponentVector = getPermutedVector(multiComponentVector, perm)
                    permutedRepresentativeIndices = getPermutedVector(representativeSitesIndices, perm)
                    permutedNumberOfAllowedSpeciesBySite = getPermutedVector(numberOfAllowedSpeciesBySite, perm)

                    clusterVectorElement += self.evaluateClusterProduct(permutedMultiComponentVector,
                                                                        permutedNumberOfAllowedSpeciesBySite,
                                                                        species_an,
                                                                        permutedRepresentativeIndices) * \
                                                                            speciesValue
                    multiplicity += speciesValue

                clusterVectorElement = clusterVectorElement/float(multiplicity)
                clusterVector.append(clusterVectorElement)
                currentMultiComponentVectorIndex += 1

        return clusterVector

    def getMultiComponentVectorPermutations(self,
                                            multiComponentVectors: list,
                                            orbitIndex: int):
        """
        @details This method return the multi-component vector permutations for each
        multi-component vector.

        Example 1: Given multi-component vectors [0, 0], [0, 1] and [1, 1]
        the returned permutations should be [[1, 0]], [[0, 1],[1, 0]], [1, 1].
        i.e. the [0, 1] multi-component vector should count elements with
        permutations [1, 0] and [1, 0].

        Example 2: Given multi-component vectors [0, 0], [0, 1], [1, 0] and [1, 1]
        the returned permutations will only be the self permutations since the
        multi-component vectors [0, 1] and [1, 0] will handle the AB vs BA choice.

        @param multiComponentVectors multi-component vectors for this orbit
        @param orbitIndex index from which to take the allowed permutations

        @returns a vector of a vector of a vector of ints; here the innermost index
        """
        allowedPermutations = self._orbitList.orbits[orbitIndex].get_allowed_sites_permutations_vec()

        elementPermutations = []
        selfPermutations = []

        selfPermutations = copy.deepcopy(multiComponentVectors[0])

        for mc in multiComponentVectors:
            mcPermutations = []
            mcPermutations.append(selfPermutations)
            takenPermutations = []
            takenPermutations.append(selfPermutations)

            for perm in allowedPermutations:
                permutedMultiComponentVector = getPermutedVector(mc, perm)
                if (permutedMultiComponentVector not in multiComponentVectors and \
                        permutedMultiComponentVector not in takenPermutations and \
                        mc != permutedMultiComponentVector):
                    mcPermutations.append(perm)
                    takenPermutations.append(permutedMultiComponentVector)
            mcPermutations.sort()
            elementPermutations.append(mcPermutations)
        return elementPermutations

    def evaluateClusterFunction(self,
                                numberOfAllowedSpecies: int,
                                clusterFunction: int,
                                species: int):
        """
        @details Evaluates the cluster function using the specified parameters.

        The cluster functions (also "orthogonal point functions") are defined as

        .. math::

           \Theta_{n}(\sigma_p) = \begin{cases}
              1                                    &\quad \text{if}~n=0 \\
              -\cos\left(\pi(n+1)\sigma_p/M\right) &\quad \text{if n is odd} \\
              -\sin\left(\pi n   \sigma_p/M\right) &\quad \text{if n is even}
            \end{cases}

        @param numberOfAllowedSpecies number of allowed species on the site in question
        @param clusterFunction index of cluster function
        @param species index of species

        @returns the value of the cluster function
        """
        value = 2.0 * np.pi * float(int((clusterFunction + 2) / 2)) * float(species) / float(numberOfAllowedSpecies)
        return -np.cos(value) if (clusterFunction + 2)%2  == 0 else -np.sin(value)

    def evaluateClusterProduct(self,
                               multiComponentVector: list,
                               numberOfAllowedSpecies: list,
                               species: list,
                               indices: list):
        """
        @details Evaluates the full cluster product of the entire cluster.

        @param multiComponentVector multi-component vector, each element of the vector gives the index of a cluster function
        @param numberOfAllowedSpecies number of species allowed on the sites in this cluster (all sites involved are assumed to have the same number of allowed species)
        @param species species that occupy (decorate) the cluster identified by atomic number
        @param indices representative lattice indices of the cluster being computed

        @returns the cluster product
        """
        clusterProduct = 1
        for i in range(len(species)):
            maploc = self._speciesMap[indices[i]]
            locc = maploc[species[i]]
            clusterProduct *= self.evaluateClusterFunction(numberOfAllowedSpecies[i],
                                                           multiComponentVector[i],
                                                           locc)
        return clusterProduct

    def getNumberOfAllowedSpeciesBySite(self,
                                        structure: Structure, 
                                        latticeSites: list):
        """
        @details Returns the number of species allowed on each site of the provided structure.

        @param structure an atomic configuration
        @param latticeSites a list of sites

        @returns the number of allowed species for each site
        """
        numberOfAllowedSpecies = \
            [structure.get_number_of_allowed_species_by_site(latsite.index)
                for latsite in latticeSites]
        return numberOfAllowedSpecies

    def precomputeMultiComponentVectors(self):
        """
        Precomputes permutations and multicomponent vectors of each orbit
        """

        self._clusterSpaceInfo = []
        emptyVec = [0]

        self._clusterSpaceInfo.append((-1, emptyVec))
        self._multiComponentVectors = [None]*len(self._orbitList)
        self._sitePermutations = [None] *len(self._orbitList)

        for i_orbit in range(len(self)):
            permutedMCVector = []
            numberOfAllowedSpecies = self.getNumberOfAllowedSpeciesBySite(
                self._primitiveStructure, self._orbitList.orbits[i_orbit].get_representative_sites())

            multiComponentVectors = self._orbitList.orbits[i_orbit].get_mc_vectors(numberOfAllowedSpecies)
            if all([k==2 for k in numberOfAllowedSpecies]):
                sitePermutations = self.getMultiComponentVectorPermutations(multiComponentVectors, i_orbit)
                self._sitePermutations[i_orbit] = sitePermutations
                self._multiComponentVectors[i_orbit] = multiComponentVectors

            self._clusterSpaceInfo.append([(i_orbit, mcv) for mcv in multiComponentVectors])

    def getClusterSpaceInfo(self, index: int):
        """
        @details Returns a pair where pair.first is the  index of the underlying orbit in _orbitList.
                  and pair.second is the multicomponent vector for the "outerlying" orbit with index `index`
        @param index index for the "outerlying" orbit

        @todo think about a better word than outerlying.
        """
        if index >= len(self._clusterSpaceInfo):
            msg = "Out of range in ClusterSpace.getClusterSpaceInfo: "
            msg += str(index) + " >= " + str(len(self._clusterSpaceInfo))
            raise RuntimeError(msg)

        return self._clusterSpaceInfo[index];

    def pruneOrbitList(self, indices: list):
        """
        @details This function removes orbits from the underlying orbit list.
        @param indices list of orbit indices 
        """
        indices.sort()
        for i in range(len(indices), -1, -1):
            self._orbitList.removeOrbit(indices[i])

        self.precomputeMultiComponentVectors()

def getPermutedVector(v: list, indices: list):
    """
    Return the permutation of v using the permutation in indices
    """
    if len(v) != len(indices):
        raise RuntimeError("Error: vectors are not of the same size in function getPermutedVector in Symmetry.hpp")

    return [v[i] for i in indices]
