#include "ClusterSpace.hpp"

/**
@details Initialize a ClusterSpace object.
@param numberOfAllowedComponents number of allowed components for each site of the primitive structure
@param chemicalSymbols chemical symbol for each site
@param primitiveOrbitList list of orbits for the primitive structure
*/
ClusterSpace::ClusterSpace(std::vector<int> numberOfAllowedComponents,
                           std::vector<std::string> chemicalSymbols,
                           const OrbitList primitiveOrbitList)
{
    _numberOfAllowedComponents = numberOfAllowedComponents;
    _primitiveOrbitList = primitiveOrbitList;
    _primitiveStructure = primitiveOrbitList.getPrimitiveStructure();
    _primitiveStructure.setNumberOfAllowedComponents(_numberOfAllowedComponents);
    setupElementMap(chemicalSymbols);
    _isClusterSpaceInitialized = false;
};

/**
@brief Sets up a map between chemical elements and the internal species enumeration scheme.
@param elements list of chemical symbols
*/
void ClusterSpace::setupElementMap(const std::vector<std::string> elements)
{
    // convert chemical symbols to atomic numbers
    std::vector<int> intElements;
    for (const auto el : elements)
    {
        intElements.push_back(PeriodicTable::strInt[el]);
    }

    // sort list of elements
    std::sort(intElements.begin(), intElements.end());
    _elements = intElements;

    // map atomic numbers to internal species enumeration scheme
    for (size_t i = 0; i < elements.size(); i++)
    {
        _elementMap[intElements[i]] = i;
    }
}



/**
@details Calculates the cluster vector for the input structure in the cluster space.
The first element in the cluster vector will always be one (1) corresponding to the zerolet.
The remaining elements of the cluster vector represent averages over orbits
(symmetry equivalent clusters) of increasing order and size.

The length of the cluster vector is `1 + sum_i(orbit_i * orbit_i.multiComponentVectors.size())`
@param structure input configuration
@todo review the necessity to have the orderIntact argument to e.g., countOrbitList.
**/

std::vector<double> ClusterSpace::getClusterVector(const Structure &structure) const
{
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_primitiveOrbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();

    // create local orbit list and count associated clusters
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto local_orbit_list = localOrbitListGenerator.getLocalOrbitList(i);
        clusterCounts.countOrbitList(structure, local_orbit_list, orderIntact);
    }

    // check that the number of unique offsets equals the number of unit cells the supercell is made of
    int numberOfUnitcellRepetitions = structure.size() / _primitiveStructure.size();
    if (uniqueOffsets != numberOfUnitcellRepetitions)
    {
        std::string errorMessage = "The number of unique offsets does not match the number of primitive units in the input structure";
        errorMessage += " {" + std::to_string(uniqueOffsets) + "}";
        errorMessage += " != ";
        errorMessage += " {" + std::to_string(numberOfUnitcellRepetitions) + "}";
        throw std::runtime_error(errorMessage);
    }

    //
    const auto clusterMap = clusterCounts.getClusterCounts();
    std::vector<double> clusterVector;
    clusterVector.push_back(1);

// CONTINUE HERE

    // Finally begin occupying the cluster vector
    for (size_t i = 0; i < _primitiveOrbitList.size(); i++)
    {
        auto repCluster = _primitiveOrbitList.getOrbit(i).getRepresentativeCluster();
        std::vector<int> allowedOccupations;
        try {
                allowedOccupations = getNumberOfAllowedComponentsBySite(_primitiveStructure, _primitiveOrbitList.getOrbit(i).getRepresentativeSites());
            }
        catch (const std::exception& e)
        {
            throw std::runtime_error("Failed getting allowed occupations in genereteClusterVector");
        }

        // Skip rest if any sites aren't active sites (i.e. allowed occupation < 2)
        if (std::any_of(allowedOccupations.begin(), allowedOccupations.end(),[](int allowedOccupation ){ return allowedOccupation < 2; }))
        {
            continue;
        }

        auto mcVectors = _primitiveOrbitList.getOrbit(i).getMCVectors(allowedOccupations);
        auto allowedPermutationsSet = _primitiveOrbitList.getOrbit(i).getAllowedSitesPermutations();
        auto elementPermutations = getMultiComponentVectorPermutations(mcVectors, i);
        repCluster.setClusterTag(i);
        int currentMCVectorIndex = 0;
        for (const auto &mcVector : mcVectors)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            for (const auto &elementsCountPair : clusterMap.at(repCluster))
            {

                // TODO check if allowedOccupations should be permuted as well.
                for (const auto &perm : elementPermutations[currentMCVectorIndex])
                {
                    auto permutedMCVector = icet::getPermutedVector(mcVector, perm);
                    auto permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);
                    clusterVectorElement += getClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            clusterVectorElement /= ((double)multiplicity);
            clusterVector.push_back(clusterVectorElement);

            currentMCVectorIndex++;
        }
    }
    return clusterVector;
}

/// Returns the native clusters count in this structure, i.e. only clusters inside the unit cell
ClusterCounts ClusterSpace::getNativeClusters(const Structure &structure) const
{
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_primitiveOrbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();
    int tags = 0;
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto local_orbit_list = localOrbitListGenerator.getLocalOrbitList(i);
        for (int j = 0; j < local_orbit_list.size(); j++)
        {
            Cluster repr_cluster = local_orbit_list.getOrbit(j).getRepresentativeCluster();

            for (const auto sites : local_orbit_list.getOrbit(j).getPermutedEquivalentSites())
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
                repr_cluster.setClusterTag(j);
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

  @param mcVectors the mc vectors for this orbit
  @param orbitIndex : The orbit index to take the allowed permutations from.

*/

std::vector<std::vector<std::vector<int>>> ClusterSpace::getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &mcVectors, const int orbitIndex) const
{
    const auto allowedPermutations = _primitiveOrbitList.getOrbit(orbitIndex).getAllowedSitesPermutations();

    std::vector<std::vector<std::vector<int>>> elementPermutations;
    std::vector<int> selfPermutation;
    for (int i = 0; i < mcVectors[0].size(); i++)
    {
        selfPermutation.push_back(i);
    }

    for (const auto &mc : mcVectors)
    {
        std::vector<std::vector<int>> mcPermutations;
        mcPermutations.push_back(selfPermutation);
        std::vector<std::vector<int>> takenPermutations;
        takenPermutations.push_back(selfPermutation);
        for (const std::vector<int> perm : allowedPermutations)
        {
            auto permutedMcVector = icet::getPermutedVector(mc, perm);
            auto findPerm = std::find(mcVectors.begin(), mcVectors.end(), permutedMcVector);
            auto findIfTaken = std::find(takenPermutations.begin(), takenPermutations.end(), permutedMcVector);
            if (findPerm == mcVectors.end() && findIfTaken == takenPermutations.end() && mc != permutedMcVector)
            {
                mcPermutations.push_back(perm);
                takenPermutations.push_back(permutedMcVector);
            }
        }
        std::sort(mcPermutations.begin(), mcPermutations.end());
        elementPermutations.push_back(mcPermutations);
    }
    return elementPermutations;
}

/**
This is the default cluster function
*/
double ClusterSpace::defaultClusterFunction(const int numberOfAllowedComponents, const int clusterFunction, const int element) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)numberOfAllowedComponents));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)numberOfAllowedComponents));
    }
}

/// Returns the full cluster product of the entire cluster (elements vector) assuming all sites have same numberOfAllowedComponents.
double ClusterSpace::getClusterProduct(const std::vector<int> &mcVector, const std::vector<int> &numberOfAllowedComponents, const std::vector<int> &elements) const
{
    double clusterProduct = 1;
    for (int i = 0; i < elements.size(); i++)
    {
        clusterProduct *= defaultClusterFunction(numberOfAllowedComponents[i], mcVector[i], _elementMap.at(elements[i]));
    }
    return clusterProduct;
}

/// Returns the allowed occupations on the sites
std::vector<int> ClusterSpace::getNumberOfAllowedComponentsBySite(const Structure &structure, const std::vector<LatticeSite> &latticeNeighbors) const
{
    std::vector<int> numberOfAllowedComponents;
    numberOfAllowedComponents.reserve(latticeNeighbors.size());
    for (const auto &latnbr : latticeNeighbors)
    {
        numberOfAllowedComponents.push_back(structure.getNumberOfAllowedComponentsBySite(latnbr.index()));
    }
    return numberOfAllowedComponents;
}

/// Collect information about the cluster space
void ClusterSpace::collectClusterSpaceInfo()
{
    _clusterSpaceInfo.clear();
    std::vector<int> emptyVec = {0};
    _clusterSpaceInfo.push_back(std::make_pair(-1, emptyVec));
    for (int i = 0; i < _primitiveOrbitList.size(); i++)
    {
        auto allowedOccupations = getNumberOfAllowedComponentsBySite(_primitiveStructure, _primitiveOrbitList.getOrbit(i).getRepresentativeSites());
        auto mcVectors = _primitiveOrbitList.getOrbit(i).getMCVectors(allowedOccupations);
        for (const auto &mcVector : mcVectors)
        {
            _clusterSpaceInfo.push_back(std::make_pair(i, mcVector));
        }
    }
    _isClusterSpaceInitialized = true;
}

/// Returns information about the cluster space
std::pair<int, std::vector<int>> ClusterSpace::getClusterSpaceInfo(const unsigned int index)
{
    if (!_isClusterSpaceInitialized)
    {
        collectClusterSpaceInfo();
    }

    if (index >= _clusterSpaceInfo.size())
    {
        std::string errMSG = "Out of range in ClusterSpace::getClusterSpaceInfo " + std::to_string(index) + " >= " + std::to_string(_clusterSpaceInfo.size());
        throw std::out_of_range(errMSG);
    }

    return _clusterSpaceInfo[index];
}
