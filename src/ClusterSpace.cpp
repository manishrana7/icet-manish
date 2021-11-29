#define _USE_MATH_DEFINES
#include <cmath>

#include "ClusterSpace.hpp"

/**
@details This constructor initializes a ClusterSpace object.
@param chemicalSymbols vector of allowed chemical symbol for each site
@param orbitList list of orbits for the primitive structure
@param positionTolerance tolerance applied when comparing positions in Cartesian coordinates
@param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates
*/
ClusterSpace::ClusterSpace(std::shared_ptr<OrbitList> orbitList,
                           const double positionTolerance,
                           const double fractionalPositionTolerance)
    : _primitiveOrbitList(orbitList)
{
    // Set up a map between atomic numbers and the internal species enumeration scheme.
    for (const auto &atomicNumbers : _primitiveOrbitList->structure().allowedAtomicNumbers())
    {
        std::unordered_map<int, int> speciesMap;
        std::vector<int> atomicNumbersCopy = atomicNumbers;
        sort(atomicNumbersCopy.begin(), atomicNumbersCopy.end());
        for (size_t i = 0; i < atomicNumbersCopy.size(); i++)
        {
            speciesMap[atomicNumbersCopy[i]] = i;
        }
        _speciesMaps.push_back(speciesMap);
    }
}

/**
@details This function calculates and then returns the cluster vector for the
input structure in the cluster space.

The first element in the cluster vector will always be one (1) corresponding to
the zerolet. The remaining elements of the cluster vector represent averages
over orbits (symmetry equivalent clusters) of increasing order and size.

@param structure input configuration
@param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates

@todo review the necessity for having the keepOrder argument to e.g., countOrbitList.
**/

std::vector<double> ClusterSpace::getClusterVector(const Structure &structure,
                                                   const double fractionalPositionTolerance) const
{
    // Construct orbit list for this structure.
    std::shared_ptr<Structure> supercell = std::make_shared<Structure>(structure);
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(*_primitiveOrbitList, supercell, fractionalPositionTolerance);
    auto currentOrbitList = localOrbitListGenerator.getFullOrbitList();

    // Check that the number of unique offsets equals the number of unit cells in the supercell.
    if (localOrbitListGenerator.getNumberOfUniqueOffsets() != structure.size() / _primitiveOrbitList->structure().size())
    {
        std::ostringstream msg;
        msg << "The number of unique offsets does not match the number of primitive units in the input structure (ClusterSpace::getClusterVector)" << std::endl;
        msg << localOrbitListGenerator.getNumberOfUniqueOffsets() << " != " << structure.size() / _primitiveOrbitList->structure().size();
        throw std::runtime_error(msg.str());
    }
    return getClusterVectorFromOrbitList(currentOrbitList, supercell);
}

/**
  @details This method returns the multi-component vector permutations for
  each multi-component vector.

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

  @todo This function should take an Orbit rather than an orbit index.
*/

/*
std::vector<std::vector<std::vector<int>>> ClusterSpace::getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &multiComponentVectors, const int orbitIndex) const
{
    const auto allowedPermutations = _primitiveOrbitList->getOrbit(orbitIndex).getAllowedClusterPermutations();

    std::vector<std::vector<std::vector<int>>> elementPermutations;
    std::vector<int> selfPermutation;
    for (size_t i = 0; i < multiComponentVectors[0].size(); i++)
    {
        selfPermutation.push_back(i);
    }

    for (const auto &mc : multiComponentVectors)
    {
        std::vector<std::vector<int>> mcPermutations;
        mcPermutations.push_back(selfPermutation);
        std::vector<std::vector<int>> takenPermutations;
        takenPermutations.push_back(selfPermutation);
        for (const std::vector<int> perm : allowedPermutations)
        {
            auto permutedMultiComponentVector = icet::getPermutedVector(mc, perm);
            auto findPerm = find(multiComponentVectors.begin(), multiComponentVectors.end(), permutedMultiComponentVector);
            auto findIfTaken = find(takenPermutations.begin(), takenPermutations.end(), permutedMultiComponentVector);
            if (findPerm == multiComponentVectors.end() && findIfTaken == takenPermutations.end() && mc != permutedMultiComponentVector)
            {
                mcPermutations.push_back(perm);
                takenPermutations.push_back(permutedMultiComponentVector);
            }
        }
        sort(mcPermutations.begin(), mcPermutations.end());
        elementPermutations.push_back(mcPermutations);
    }
    return elementPermutations;
}
*/

/**
@details Evaluates the cluster function using the specified parameters.

The cluster functions (also "orthogonal point functions") are defined as

.. math::

   \\Theta_{n}(\\sigma_p) = \\begin{cases}
      1                                     &\\quad \\text{if}~n=0 \\
      -\\cos\\left(\\pi(n+1)\\sigma_p/M\\right) &\\quad \\text{if n is odd} \\
      -\\sin\\left(\\pi n   \\sigma_p/M\\right) &\\quad \\text{if n is even}
    \\end{cases}

@param numberOfAllowedSpecies number of allowed species on the site in question
@param clusterFunction index of cluster function
@param species index of species

@returns the value of the cluster function
*/
double ClusterSpace::evaluateClusterFunction(const int numberOfAllowedSpecies, const int clusterFunction, const int species) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)species / ((double)numberOfAllowedSpecies));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)species / ((double)numberOfAllowedSpecies));
    }
}

/**
@details Evaluates the full cluster product of the entire cluster.

@param multiComponentVector multi-component vector, each element of the vector gives the index of a cluster function
@param numberOfAllowedSpecies number of species allowed on the sites in this cluster (all sites involved are assumed to have the same number of allowed species)
@param species species that occupy (decorate) the cluster identified by atomic number
@param indices representative lattice indices of the cluster being computed
@param permutation Describes the desired order of the points in the cluster; for example,
        if the current cluster is a triplet cluster and permutation = {1, 0, 2}, then multiComponentVector, numberOfAllowedSpecies, and indices should be permuted
        such that their first two elements change places

@returns the cluster product
**/
double ClusterSpace::evaluateClusterProduct(const std::vector<int> &multiComponentVector, const std::vector<int> &numberOfAllowedSpecies, const std::vector<int> &species, const std::vector<int> &indices, const std::vector<int> &permutation) const
{
    double clusterProduct = 1;

    for (size_t i = 0; i < species.size(); i++)
    {
        int index = permutation[i];
        clusterProduct *= evaluateClusterFunction(numberOfAllowedSpecies[index], multiComponentVector[index], _speciesMaps[indices[index]].at(species[i]));
    }
    return clusterProduct;
}

/// Computes permutations and multicomponent vectors of each orbit.
/*
void ClusterSpace::computeMultiComponentVectors()
{
    std::vector<int> emptyVec = {0};
    _clusterVectorElementInfoList.clear();
    _clusterVectorElementInfoList.resize(_primitiveOrbitList->size());

    int clusterVectorIndex = 0;
    for (size_t orbitIndex = 0; orbitIndex < _primitiveOrbitList->size(); orbitIndex++)
    {

        auto numberOfAllowedSpecies = _primitiveStructure->getNumberOfAllowedSpeciesBySites(_primitiveOrbitList->getOrbit(orbitIndex).representativeCluster().latticeSites());
        auto multiComponentVectors = _primitiveOrbitList->getOrbit(orbitIndex).getMultiComponentVectors(numberOfAllowedSpecies);
        if (std::none_of(numberOfAllowedSpecies.begin(), numberOfAllowedSpecies.end(), [](int n)
                         { return n < 2; }))
        {
            auto sitePermutations = getMultiComponentVectorPermutations(multiComponentVectors, orbitIndex);
            for (int j = 0; j < multiComponentVectors.size(); j++)
            {
                clusterVectorIndex++;
                double multiplicity = (double)sitePermutations[j].size() * (double)_primitiveOrbitList->getOrbit(orbitIndex).size() / (double)_primitiveStructure->size();
                ClusterVectorElementInfo cvInfo = {multiComponentVectors[j],
                                                   sitePermutations[j],
                                                   clusterVectorIndex,
                                                   multiplicity};
                _clusterVectorElementInfoList[orbitIndex].push_back(cvInfo);
            }
        }
    }
    _clusterVectorLength = clusterVectorIndex + 1;
}
*/

size_t ClusterSpace::size() const
{
    int size = 1; // 1 for the zerolet
    for (size_t orbitIndex = 0; orbitIndex < _primitiveOrbitList->size(); orbitIndex++)
    {
        const Orbit &orbit = _primitiveOrbitList->getOrbit(orbitIndex);
        size += orbit._clusterVectorElements.size();
    }
    return size;
}

/**
@details This function removes orbits from the underlying orbit list.
@param indices list of orbit indices
**/
void ClusterSpace::removeOrbits(std::vector<size_t> &indices)
{
    std::sort(indices.begin(), indices.end());
    for (int i = indices.size() - 1; i >= 0; i--)
    {
        _primitiveOrbitList->removeOrbit(indices[i]);
    }
}

/*
@details Occupy cluster vector based on a supercell and a corresponding orbit list.
@param orbitList An orbit list to be used for counting. Can be either a full orbit list or a "loca"
 orbit list
@param supercell Defines the occupations of the structure whose cluster vector should be computed
@param firstElement First element of the cluster vector (default: 1.0)
@param flipIndex If a local cluster vector should be calculated this argument is used to specify the index of the site whose local cluster vector should be computed. If calculating a change in cluster vector, this is the site whose occupation has changed. If -1 (default), the total cluster vector will be calculated.
@param newOccupation New atomic number on the site with index flipIndex. If this argument is not -1, a change in cluster vector will be calculated.
*/
const std::vector<double> ClusterSpace::getClusterVectorFromOrbitList(const OrbitList &orbitList,
                                                                      const std::shared_ptr<Structure> supercell,
                                                                      const double firstElement,
                                                                      const int flipIndex,
                                                                      const int newOccupation) const
{
    if (newOccupation >= 0 && flipIndex == -1)
    {
        throw std::runtime_error("flipIndex needs to be specified (larger than -1) if newOccupation is specified (ClusterSpace::getClusterVectorFromOrbitList)");
    }
    std::vector<double> clusterVector;
    clusterVector.push_back(firstElement);
    if (_primitiveOrbitList->size() != orbitList.size())
    {
        std::cout << orbitList.size() << " >= " << _primitiveOrbitList->size() << std::endl;
        throw std::runtime_error("Orbit lists do not match (ClusterSpace::getClusterVectorFromOrbitList)");
    }

    for (size_t currentOrbitIndex = 0; currentOrbitIndex < _primitiveOrbitList->size(); currentOrbitIndex++)
    {
        const Orbit &currentOrbit = orbitList.getOrbit(currentOrbitIndex);
        const Orbit &currentPrimitiveOrbit = _primitiveOrbitList->getOrbit(currentOrbitIndex);

        // Count clusters
        std::map<std::vector<int>, double> counts;
        if (newOccupation > -1)
        {
            counts = currentOrbit.getClusterCountChanges(supercell, flipIndex, newOccupation);
        }
        else
        {
            counts = currentOrbit.getClusterCounts(supercell, flipIndex);
        }

        // Extract allowed occupations
        // @todo Get this from Orbit instead
        std::vector<int> allowedOccupations = _primitiveOrbitList->structure().getNumberOfAllowedSpeciesBySites(currentPrimitiveOrbit.representativeCluster().latticeSites());

        // Skip the rest if any of the sites are inactive (i.e. allowed occupation < 2)
        // @todo Let orbit handle this
        if (std::any_of(allowedOccupations.begin(), allowedOccupations.end(), [](int allowedOccupation)
                        { return allowedOccupation < 2; }))
        {
            continue;
        }

        std::vector<int> indicesOfRepresentativeSites;
        for (const LatticeSite &site : currentPrimitiveOrbit.representativeCluster().latticeSites())
        {
            indicesOfRepresentativeSites.push_back(site.index());
        }

        // Loop over all multi-component vectors for this orbit.
        // These are vectors of ints (where the int represents a cluster function index).
        // Example 1: For an AB alloy we obtain [0, 0] and [0, 0, 0] for pair and triplet terms, respectively.
        // Example 2: For an ABC alloy we obtain [0, 0], [0, 1], [1, 1] for pairs and similarly for triplets.
        // Depending on the symmetry of the cluster one might also obtain [1, 0] (e.g., in a clathrate or for some clusters on a HCP lattice).
        for (auto &cvElement : currentPrimitiveOrbit._clusterVectorElements)
        {
            /// Push back zero if nothing was counted for this orbit
            if (counts.size() == 0)
            {
                clusterVector.push_back(0);
                continue;
            }

            double clusterVectorElement = 0;

            /// Loop over all the counts for this orbit
            for (const auto &elementsCountPair : counts)
            {
                /// Loop over all equivalent permutations for this orbit and mc vector
                for (const auto &perm : cvElement.sitePermutations)
                {
                    clusterVectorElement += evaluateClusterProduct(cvElement.multicomponentVector,
                                                                   allowedOccupations,
                                                                   elementsCountPair.first,
                                                                   indicesOfRepresentativeSites,
                                                                   perm) *
                                            elementsCountPair.second;
                }
            }

            // Usually we could have counted multiplicity by simply adding the number of
            // clusters in the orbit (elementsCountPair.second), but in the case of
            // local cluster vectors or changes in cluster vectors, we have only counted
            // a subset of the clusters. We thus need to compute the multiplicity by
            // analyzing the orbit in detail.
            clusterVector.push_back(clusterVectorElement / (double)cvElement.multiplicity * (double)_primitiveOrbitList->structure().size() / (double)supercell->size());
        }
    }
    return clusterVector;
}
