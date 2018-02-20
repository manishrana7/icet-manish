#include "ClusterSpace.hpp"

/**
Calculate a cluster vector from the structure and the internal state of the cluster space.

The first element in the cluster vector will always be a constant (1).
The second element in the cluster vector is the average of the first orbit.
The third element might be the first orbits second multi-component vector
or it is the second orbit if it only has one multi-component vector.

The length of the cluster vector will be `1 + sum_i(orbit_i * orbit_i.multicimponentvectors.size())`
**/

std::vector<double> ClusterSpace::generateClusterVector(const Structure &structure2) const
{
    Structure structure = structure2;
    //structure.setNumberOfAllowedComponents(_Mi);
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_primitive_orbit_list, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getUniqueOffsetsCount();
    ClusterCounts clusterCounts = ClusterCounts();
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto local_orbit_list = localOrbitListGenerator.generateLocalOrbitList(i);
        clusterCounts.countOrbitList(structure, local_orbit_list, orderIntact);
    }

    if (uniqueOffsets != structure2.size() / _primitive_structure.size())
    {
        std::string errorMessage = "The number of unique offsets do not match supercell.size() / primitive.size()";
        errorMessage += " {" + std::to_string(uniqueOffsets) + "}";
        errorMessage += " != ";
        errorMessage += " {" + std::to_string(structure2.size() / _primitive_structure.size()) + "}";
        throw std::runtime_error(errorMessage);
    }

    const auto clusterMap = clusterCounts.getClusterCounts();
    std::vector<double> clusterVector;
    clusterVector.push_back(1);
    // Finally begin occupying the cluster vector
    for (size_t i = 0; i < _primitive_orbit_list.size(); i++)
    {
        auto repCluster = _primitive_orbit_list.getOrbit(i).getRepresentativeCluster();
        auto allowedOccupations = getAllowedOccupations(_primitive_structure, _primitive_orbit_list.getOrbit(i).getRepresentativeSites());
        auto mcVectors = _primitive_orbit_list.getOrbit(i).getMCVectors(allowedOccupations);
        auto allowedPermutationsSet = _primitive_orbit_list.getOrbit(i).getAllowedSitesPermutations();
        auto elementPermutations = getMCVectorPermutations(mcVectors, i);
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
                    auto permutedMCVector = icet::getPermutatedVector(mcVector, perm);
                    clusterVectorElement += getClusterProduct(permutedMCVector, allowedOccupations, elementsCountPair.first) * elementsCountPair.second;
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
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_primitive_orbit_list, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getUniqueOffsetsCount();
    ClusterCounts clusterCounts = ClusterCounts();
    int tags = 0;
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto local_orbit_list = localOrbitListGenerator.generateLocalOrbitList(i);
        for (int j = 0; j < local_orbit_list.size(); j++)
        {
            Cluster repr_cluster = local_orbit_list.getOrbit(j).getRepresentativeCluster();

            for (const auto sites : local_orbit_list.getOrbit(j).getPermutatedEquivalentSites())
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

std::vector<std::vector<std::vector<int>>> ClusterSpace::getMCVectorPermutations(const std::vector<std::vector<int>> &mcVectors, const int orbitIndex) const
{
    const auto allowedPermutations = _primitive_orbit_list.getOrbit(orbitIndex).getAllowedSitesPermutations();

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
            auto permutedMcVector = icet::getPermutatedVector(mc, perm);
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
double ClusterSpace::defaultClusterFunction(const int Mi, const int clusterFunction, const int element) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)Mi));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)Mi));
    }
}

/// Returns the full cluster product of the entire cluster (elements vector) assuming all sites have same Mi.
double ClusterSpace::getClusterProduct(const std::vector<int> &mcVector, const std::vector<int> &Mi, const std::vector<int> &elements) const
{
    double clusterProduct = 1;
    for (int i = 0; i < elements.size(); i++)
    {
        clusterProduct *= defaultClusterFunction(Mi[i], mcVector[i], _elementRepresentation.at(elements[i]));
    }
    return clusterProduct;
}

/// Returns the allowed occupations on the sites
std::vector<int> ClusterSpace::getAllowedOccupations(const Structure &structure, const std::vector<LatticeSite> &latticeNeighbors) const
{
    std::vector<int> Mi;
    Mi.reserve(latticeNeighbors.size());
    for (const auto &latnbr : latticeNeighbors)
    {
        Mi.push_back(structure.getNumberOfAllowedComponents(latnbr.index()));
    }
    return Mi;
}

/// Collect information about the cluster space
void ClusterSpace::setupClusterSpaceInfo()
{
    _clusterSpaceInfo.clear();
    std::vector<int> emptyVec = {0};
    _clusterSpaceInfo.push_back(std::make_pair(-1, emptyVec));
    for (int i = 0; i < _primitive_orbit_list.size(); i++)
    {
        auto allowedOccupations = getAllowedOccupations(_primitive_structure, _primitive_orbit_list.getOrbit(i).getRepresentativeSites());
        auto mcVectors = _primitive_orbit_list.getOrbit(i).getMCVectors(allowedOccupations);
        for (const auto &mcVector : mcVectors)
        {
            _clusterSpaceInfo.push_back(std::make_pair(i, mcVector));
        }
    }
    _isClusterSpaceInitialized = true;
}

/// Retrieve information about the cluster space
std::pair<int, std::vector<int>> ClusterSpace::getClusterSpaceInfo(const unsigned int index)
{
    if (!_isClusterSpaceInitialized)
    {
        setupClusterSpaceInfo();
    }

    if (index >= _clusterSpaceInfo.size())
    {
        std::string errMSG = "Out of range in ClusterSpace::getClusterSpaceInfo " + std::to_string(index) + " >= " + std::to_string(_clusterSpaceInfo.size());
        throw std::out_of_range(errMSG);
    }

    return _clusterSpaceInfo[index];
}

/// Returns the cluster space size i.e. the length of a cluster vector
size_t ClusterSpace::getClusterSpaceSize()
{
    if (!_isClusterSpaceInitialized)
    {
        setupClusterSpaceInfo();
    }
    // Plus one for zerolet
    return _clusterSpaceInfo.size();
}
