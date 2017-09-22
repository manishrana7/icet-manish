#include "ClusterCounts.hpp"

/**
Counts clusters given this compact form of latticeneighbors (see ManybodyNeighborlist for more details)
*/
// build(const Neighborlist &nl, int index, int order, bool);
void ClusterCounts::countLatticeNeighbors(const Structure &structure,
                                          const std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &latticeNeighbors)
{
    for (const auto &neighborPair : latticeNeighbors)
    {
        //Now we have std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>
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
                          const std::vector<LatticeNeighbor> &latticeNeighbors)
{
    size_t clusterSize = latticeNeighbors.size();
    std::vector<int> elements(clusterSize);
    for (size_t i = 0; i < latticeNeighbors.size(); i++)
    {
        elements[i] = structure.getElement(latticeNeighbors[i].index);
    }

    Cluster cluster = Cluster(structure, latticeNeighbors);
    countCluster(cluster, elements);
}

/**
Will count the vectors in latticeNeighbors and assuming these sets of sites are represented by the cluster 'cluster'
*/
void ClusterCounts::count(const Structure &structure, const std::vector<std::vector<LatticeNeighbor>> &latticeNeighbors,
                          const Cluster &cluster)
{

    for (const auto &latnbrs : latticeNeighbors)
    {
        std::vector<int> elements(latnbrs.size());
        for (size_t i = 0; i < latnbrs.size(); i++)
        {
            elements[i] = structure.getElement(latnbrs[i].index);
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
    Counts the clusters of the sites in each orbit of orbitlist

    Parameters
    ----------
    structure:
        icet structure object to be counted on

    orbitlist:
        icet orbitlist class
*/
void ClusterCounts::countOrbitlist(const Structure &structure, const OrbitList &orbitlist)
{
    for (int i = 0; i < orbitlist.size(); i++)
    {
        Cluster repr_cluster = orbitlist.getOrbit(i).getRepresentativeCluster();
        repr_cluster.setClusterTag(i);

        count(structure, orbitlist.getOrbit(i).getEquivalentSites(), repr_cluster);
    }
}