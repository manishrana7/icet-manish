#include "ClusterCounts.hpp"

/// Count clusters given this compact form of lattice neighbors (see ManyBodyNeighborList for more details)
// build(const NeighborList &nl, int index, int order, bool);
void ClusterCounts::countLatticeSites(const Structure &structure,
                                      const std::vector<std::pair<std::vector<LatticeSite>,std::vector<LatticeSite>>> &latticeNeighbors)
{
    for (const auto &neighborPair : latticeNeighbors)
    {
        //Now we have std::pair<std::vector<LatticeSite>, std::vector<LatticeSite>>
        //pair.first == the base indices and pair.second is all indices that form clusters with the base indices
        if (neighborPair.second.size() > 0)
        {
            for (const auto &combinationIndice : neighborPair.second)
            {
                auto latticePointsForCluster = neighborPair.first;
                latticePointsForCluster.push_back(combinationIndice);
                count(structure, latticePointsForCluster);
            }
        }
        else
        {
            //count singlets here
            count(structure, neighborPair.first);
        }
    }
}
/**
The simplest form of counting clusters using the mbnl format

Get the indice of one set of indices and counts this
*/
void ClusterCounts::count(const Structure &structure,
                          const std::vector<LatticeSite> &latticeNeighbors)
{
    size_t clusterSize = latticeNeighbors.size();
    std::vector<int> elements(clusterSize);
    for (size_t i = 0; i < latticeNeighbors.size(); i++)
    {
        elements[i] = structure.getAtomicNumber(latticeNeighbors[i].index());
    }

    // Don't do intact order since there is no reason for it
    Cluster cluster = Cluster(structure, latticeNeighbors);
    countCluster(cluster, elements, false);
}

/**
Will count the vectors in latticeNeighbors and assuming these sets of sites are represented by the cluster 'cluster'
*/
void ClusterCounts::count(const Structure &structure, const std::vector<std::vector<LatticeSite>> &latticeNeighbors,
                          const Cluster &cluster, bool orderIntact)
{

    for (const auto &latnbrs : latticeNeighbors)
    {
        std::vector<int> elements(latnbrs.size());
        for (size_t i = 0; i < latnbrs.size(); i++)
        {
            elements[i] = structure.getAtomicNumber(latnbrs[i].index());
        }
        countCluster(cluster, elements, orderIntact);
    }
}


///Count cluster only through this function
void ClusterCounts::countCluster(const Cluster &cluster, const std::vector<int> &elements, bool orderIntact)
{
    if(orderIntact)
    {
        _clusterCounts[cluster][elements] += 1;
    }
    else
    {
        std::vector<int> sortedElements = elements;
        std::sort(sortedElements.begin(), sortedElements.end());
        _clusterCounts[cluster][sortedElements ] += 1;
    }
}


/**
 @brief Counts the clusters in the input structure.
 @param structure input configuration
 @param orbitList orbit list
 @param orderIntact if True do not reorder clusters before comparison (i.e., ABC != ACB)
*/
void ClusterCounts::countOrbitList(const Structure &structure, const OrbitList &orbitList, bool orderIntact)
{
    for (int i = 0; i < orbitList.size(); i++)
    {
        Cluster repr_cluster = orbitList.getOrbit(i).getRepresentativeCluster();
        repr_cluster.setClusterTag(i);
        if(orderIntact && repr_cluster.order()!= 1)
        {
            count(structure, orbitList.getOrbit(i).getPermutedEquivalentSites(), repr_cluster, orderIntact);
        }
        else
        {
            count(structure, orbitList.getOrbit(i).getEquivalentSites(), repr_cluster, orderIntact);
        }
    }
}

// Return clusters count information (elements and count) per cluster
void ClusterCounts::setupClusterCountsInfo()
{
    _clusterCountsInfo.clear();
    for (const auto &map_pair : _clusterCounts)
    {
      for (const auto &element_count_pair : map_pair.second)
      {
        std::vector<std::string> elementVec;
        for (auto el : element_count_pair.first)
        {
            auto getElementSymbol = PeriodicTable::intStr[el];
            elementVec.push_back(getElementSymbol);
        }
        auto getCountPair = map_pair.second.at(element_count_pair.first);
        _clusterCountsInfo.push_back(std::make_pair(elementVec, getCountPair));
      }
    }
    std::sort(
        _clusterCountsInfo.begin(), 
        _clusterCountsInfo.end(),
        [](const std::pair<std::vector<std::string>, int> &a, const std::pair<std::vector<std::string>, int> &b){
            return a.first.size() < b.first.size();
        });
}

// Retrieve cluster count information
std::pair<std::vector<std::string>, int> ClusterCounts::getClusterCountsInfo(const unsigned int index)
{
    if (index >= _clusterCountsInfo.size())
    {
        std::string errMSG = "Out of range in ClusterCounts::getClusterCountsInfo " + std::to_string(index) + " >= " + std::to_string(_clusterCountsInfo.size());
        throw std::out_of_range(errMSG);
    }
    return _clusterCountsInfo[index];
}