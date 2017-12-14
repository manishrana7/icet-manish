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

    Cluster cluster = Cluster(structure, latticeNeighbors);
    countCluster(cluster, elements);
}

/**
Will count the vectors in latticeNeighbors and assuming these sets of sites are represented by the cluster 'cluster'
*/
void ClusterCounts::count(const Structure &structure, const std::vector<std::vector<LatticeSite>> &latticeNeighbors,
                          const Cluster &cluster)
{

    for (const auto &latnbrs : latticeNeighbors)
    {
        std::vector<int> elements(latnbrs.size());
        for (size_t i = 0; i < latnbrs.size(); i++)
        {
            elements[i] = structure.getAtomicNumber(latnbrs[i].index());
        }
        countCluster(cluster, elements);
    }
}

///Count cluster only through this function
void ClusterCounts::countCluster(const Cluster &cluster, const std::vector<int> &elements)
{
    std::vector<int> sortedElements = elements;
    std::sort(sortedElements.begin(), sortedElements.end());
    _clusterCounts[cluster][sortedElements] += 1;
}

/**
    Count the clusters of the sites in each orbit of the orbit list

    Parameters
    ----------
    structure:
        icet structure object to be counted on

    orbitList:
        OrbitList object
*/
void ClusterCounts::countOrbitList(const Structure &structure, const OrbitList &orbitList, bool orderIntact)
{

    for (int i = 0; i < orbitList.size(); i++)
    {
        Cluster repr_cluster = orbitList.getOrbit(i).getRepresentativeCluster();
        repr_cluster.setClusterTag(i);
        if(orderIntact && repr_cluster.order()!= 1)
        {
            count(structure, orbitList.getOrbit(i).getPermutatedEquivalentSites(), repr_cluster);
        }
        else
        {
            count(structure, orbitList.getOrbit(i).getEquivalentSites(), repr_cluster);
        }
    }

}
