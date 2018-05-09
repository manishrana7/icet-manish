#include "ClusterSpace.hpp"

//namespace icet {

/**
@details This constructor initializes a ClusterSpace object.
@param numberOfAllowedSpecies number of allowed components for each site of the primitive structure
@param chemicalSymbols chemical symbol for each site
@param orbitList list of orbits for the primitive structure
*/
ClusterSpace::ClusterSpace(std::vector<int> numberOfAllowedSpecies,
                           std::vector<std::string> chemicalSymbols,
                           const OrbitList orbitList)
{
    _numberOfAllowedSpeciesPerSite = numberOfAllowedSpecies;
    _orbitList = orbitList;
    _primitiveStructure = orbitList.getPrimitiveStructure();
    _primitiveStructure.setNumberOfAllowedSpecies(_numberOfAllowedSpeciesPerSite);

    // Set up a map between chemical elements and the internal species enumeration scheme.
    for (const auto el : chemicalSymbols)
    {
        _elements.push_back(PeriodicTable::strInt[el]);
    }
    sort(_elements.begin(), _elements.end());
    for (size_t i = 0; i < _elements.size(); i++)
    {
        _elementMap[_elements[i]] = i;
    }

    /// @todo Why is the collectClusterSpaceInfo function not executed
    /// immediately, which would render _isClusterSpaceInitialized unnecessary?
    _isClusterSpaceInitialized = false;
}


/**
@details This function calculates and then returns the cluster vector for the
input structure in the cluster space.

The first element in the cluster vector will always be one (1) corresponding to
the zerolet. The remaining elements of the cluster vector represent averages
over orbits (symmetry equivalent clusters) of increasing order and size.

@param structure input configuration

@todo review the necessity for having the orderIntact argument to e.g.,
countOrbitList.
**/
std::vector<double> ClusterSpace::getClusterVector(const Structure &structure) const
{
    // count the clusters in the orbit with the same orientation (order) as the prototype cluster
    // @todo Clarify description.
    bool orderIntact = true;

    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_orbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();

    // Create local orbit list and count associated clusters.
    for (int i = 0; i < uniqueOffsets; i++) {
        const auto localOrbitList = localOrbitListGenerator.getLocalOrbitList(i);
        clusterCounts.countOrbitList(structure, localOrbitList, orderIntact);
    }

    // Check that the number of unique offsets equals the number of unit cells in the supercell.
    int numberOfUnitcellRepetitions = structure.size() / _primitiveStructure.size();
    if (uniqueOffsets != numberOfUnitcellRepetitions) {
        std::string msg = "The number of unique offsets does not match the number of primitive units in the input structure: ";
        msg += std::to_string(uniqueOffsets) + " != " + std::to_string(numberOfUnitcellRepetitions);
        throw std::runtime_error(msg);
    }

    /// @todo Describe what is happening here.
    const auto clusterMap = clusterCounts.getClusterCounts();

    // Initialize cluster vector and insert zerolet.
    std::vector<double> clusterVector;
    clusterVector.push_back(1);

    // Loop over orbits.
    /// @todo Turn this into a proper loop over orbits. This probably requires upgrading OrbitList.
    for (size_t i = 0; i < _orbitList.size(); i++)
    {

        auto representativeCluster = _orbitList.getOrbit(i).getRepresentativeCluster();
        // @todo This is necessary. Side effects need to carefully evaluated. Ideally move this task elsewhere as it is repeated for every structure, for which the cluster vector is computed.
        representativeCluster.setTag(i);

        std::vector<int> numberOfAllowedSpeciesForEachSite;
        try {
            numberOfAllowedSpeciesForEachSite = getNumberOfAllowedSpeciesForEachSite(_primitiveStructure, _orbitList.getOrbit(i).getRepresentativeSites());
        }
        catch (const std::exception& e) {
            std::string msg = "Failed retrieving the number of allowed species in ClusterSpace::getClusterVector\n";
            msg += e.what();
            throw std::runtime_error(msg);
        }

        // Jump to the next orbit if none of the sites in the representative cluster are active (i.e. the number of allowed species on this site is less than 2).
        if (any_of(numberOfAllowedSpeciesForEachSite.begin(), numberOfAllowedSpeciesForEachSite.end(), [](int n){ return n < 2; }))
        {
            continue;
        }

        // First we obtain the multi-component vectors for this orbit, i.e. a vector
        // of vectors of int (where the int represents a cluster function index).
        // Example 1: For an AB alloy we obtain [0, 0] and [0, 0, 0] for pair and triplet terms, respectively.
        // Example 2: For an ABC alloy we obtain [0, 0], [0, 1], [1, 1] for pairs and similarly for triplets.
        // Depending on the symmetry of the cluster one might also obtain [1, 0] (e.g., in a clathrate or for some clusters on a HCP lattice).
        auto multiComponentVectors = _orbitList.getOrbit(i).getMultiComponentVectors(numberOfAllowedSpeciesForEachSite);
        // @todo Make getMultiComponentVectorPermutations take an Orbit rather than an index. Then swap the loop over int for a loop over Orbit above.
        auto elementPermutations = getMultiComponentVectorPermutations(multiComponentVectors, i);
        int currentMultiComponentVectorIndex = 0;
        for (const auto &multiComponentVector : multiComponentVectors)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            for (const auto &elementsCountPair : clusterMap.at(representativeCluster))
            {

                /// @todo Check if numberOfAllowedSpecies should be permuted as well. Is this todo still relevant?
                for (const auto &perm : elementPermutations[currentMultiComponentVectorIndex])
                {
                    auto permutedMultiComponentVector = icet::getPermutedVector(multiComponentVector, perm);
                    auto permutedNumberOfAllowedSpeciesForEachSite = icet::getPermutedVector(numberOfAllowedSpeciesForEachSite, perm);
                    clusterVectorElement += getClusterProduct(permutedMultiComponentVector, permutedNumberOfAllowedSpeciesForEachSite, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            clusterVectorElement /= ((double)multiplicity);
            clusterVector.push_back(clusterVectorElement);

            currentMultiComponentVectorIndex++;
        }
    }
    return clusterVector;
}


/// Returns the native clusters count in this structure, i.e. only clusters inside the unit cell.
ClusterCounts ClusterSpace::getNativeClusters(const Structure &structure) const
{
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_orbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();
    int tags = 0;
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto localOrbitList = localOrbitListGenerator.getLocalOrbitList(i);
        for (int j = 0; j < localOrbitList.size(); j++)
        {
            Cluster repr_cluster = localOrbitList.getOrbit(j).getRepresentativeCluster();

            for (const auto sites : localOrbitList.getOrbit(j).getPermutedEquivalentSites())
            {
                bool cont = true;
                for (auto site : sites)
                {
                    if (site.unitcellOffset().norm() > 0.1)
                    {
                        cont = false;
                    }
                }
                if (!cont)
                {
                    continue;
                }
                repr_cluster.setTag(j);
                if (repr_cluster.order() != 1)
                {
                    std::vector<int> elements(sites.size());
                    for (size_t i = 0; i < sites.size(); i++)
                    {
                        elements[i] = structure.getAtomicNumber(sites[i].index());
                    }
                    clusterCounts.countCluster(repr_cluster, elements, orderIntact);
                }
                else
                {
                    std::vector<int> elements(sites.size());
                    for (size_t i = 0; i < sites.size(); i++)
                    {
                        elements[i] = structure.getAtomicNumber(sites[i].index());
                    }
                    clusterCounts.countCluster(repr_cluster, elements, orderIntact);
                }
            }
        }
    }
    return clusterCounts;
}

/**
  @details This method return the mc vector permutations for each mc vector.
  Example1: Given mc vectors [0, 0], [0,1] and [1,1]
  the returned permutations should be [[1,0]], [[0,1],[1,0]], [1,1].
  i.e. the [0,1] mc vector should count elements with permutations [1,0] and [1,0]

  Given mc vectors [0, 0], [0,1], [1,0] and [1,1] the returned permutations
  will only be the self permutations since the mc vectors [0,1] and [1,0] will handle
  the AB vs BA choice.

  @param multiComponentVectors the mc vectors for this orbit
getNumberOfAllowedSpeciesForEachSitet index to take the allowed permutations from.
getNumberOfAllowedSpeciesForEachSite

  @todo This function should take an Orbit rather than an orbit index.
  @returns a vector of a vector of a vector of ints; here the innermost index
*/

std::vector<std::vector<std::vector<int>>> ClusterSpace::getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &multiComponentVectors, const int orbitIndex) const
{
    const auto allowedPermutations = _orbitList.getOrbit(orbitIndex).getAllowedSitesPermutations();

    std::vector<std::vector<std::vector<int>>> elementPermutations;
    std::vector<int> selfPermutation;
    for (int i = 0; i < multiComponentVectors[0].size(); i++)
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

/**
This is the default cluster function
*/
double ClusterSpace::defaultClusterFunction(const int numberOfAllowedSpecies, const int clusterFunction, const int element) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)numberOfAllowedSpecies));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)numberOfAllowedSpecies));
    }
}

/// Returns the full cluster product of the entire cluster (elements vector) assuming all sites have same numberOfAllowedSpecies.
double ClusterSpace::getClusterProduct(const std::vector<int> &MultiComponentVector, const std::vector<int> &numberOfAllowedSpecies, const std::vector<int> &elements) const
{
    double clusterProduct = 1;
    for (int i = 0; i < elements.size(); i++)
    {
        clusterProduct *= defaultClusterFunction(numberOfAllowedSpecies[i], MultiComponentVector[i], _elementMap.at(elements[i]));
    }
    return clusterProduct;
}

/// Returns the number of allowed species on the sites
std::vector<int> ClusterSpace::getNumberOfAllowedSpeciesForEachSite(const Structure &structure, const std::vector<LatticeSite> &latticeNeighbors) const
{
    std::vector<int> numberOfAllowedSpecies;
    numberOfAllowedSpecies.reserve(latticeNeighbors.size());
    for (const auto &latnbr : latticeNeighbors)
    {
        numberOfAllowedSpecies.push_back(structure.getNumberOfAllowedSpeciesForEachSite(latnbr.index()));
    }
    return numberOfAllowedSpecies;
}

/// Collect information about the cluster space.
void ClusterSpace::collectClusterSpaceInfo()
{
    _clusterSpaceInfo.clear();
    std::vector<int> emptyVec = {0};
    _clusterSpaceInfo.push_back(make_pair(-1, emptyVec));
    for (int i = 0; i < _orbitList.size(); i++)
    {
        auto numberOfAllowedSpecies = getNumberOfAllowedSpeciesForEachSite(_primitiveStructure, _orbitList.getOrbit(i).getRepresentativeSites());
        auto multiComponentVectors = _orbitList.getOrbit(i).getMultiComponentVectors(numberOfAllowedSpecies);
        for (const auto &MultiComponentVector : multiComponentVectors)
        {
            _clusterSpaceInfo.push_back(make_pair(i, MultiComponentVector));
        }
    }
    _isClusterSpaceInitialized = true;
}

/// Returns information about the cluster space.
std::pair<int, std::vector<int>> ClusterSpace::getClusterSpaceInfo(const unsigned int index)
{
    if (!_isClusterSpaceInitialized)
    {
        collectClusterSpaceInfo();
    }

    if (index >= _clusterSpaceInfo.size())
    {
        std::string msg = "Out of range in ClusterSpace::getClusterSpaceInfo: ";
        msg += std::to_string(index) + " >= " + std::to_string(_clusterSpaceInfo.size());
        throw std::out_of_range(msg);
    }

    return _clusterSpaceInfo[index];
}

//}