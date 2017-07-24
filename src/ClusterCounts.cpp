#include "ClusterCounts.hpp"

///
// void ClusterCounts::count_using_mbnl(const Structure &structure, ManybodyNeighborlist &mbnl, const int order)
// {
//     bool bothways = false;

//     for (size_t latticeIndex = 0; latticeIndex < structure.size(); latticeIndex++)
//     {
//         auto latticeNeighbors = mbnl.build(structure, latticeIndex, order, bothways);
//         void countLatticeNeighbors(structure, latticeNeighbors);
//     }
// }

/**
Counts all the pairs using the neigbhorlist

*/

void ClusterCounts::count_pairs(const Structure &structure, const Neighborlist &neighborlist)
{   
    Vector3d zeroVector = {0.0, 0.0, 0.0};

    std::vector<LatticeNeighbor> pairNeighbor(2);
    pairNeighbor[0] = LatticeNeighbor(0,zeroVector);
    for(size_t i=0; i < structure.size(); i++)
    {
        auto i_neighbors = neighborlist.getNeighbors(i);
        pairNeighbor[0].index = i;
        for(const auto &neighbor : i_neighbors)
        {
          pairNeighbor[1] = neighbor;
          count(structure, pairNeighbor);
        }

    }

}


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
        for (const auto &combinationIndice : neighborPair.second)
        {
            auto latticePointsForCluster = neighborPair.first;
            latticePointsForCluster.push_back(combinationIndice);
            count(structure, latticePointsForCluster);
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
    std::vector<int> sites(clusterSize);
    std::vector<double> distances;
    std::vector<int> elements(clusterSize);
    distances.reserve((clusterSize * (clusterSize - 1) / 2));
    for (size_t i = 0; i < latticeNeighbors.size(); i++)
    {
        sites[i] = structure.getSite(latticeNeighbors[i].index);
        elements[i] = structure.getElement(latticeNeighbors[i].index);
        for (size_t j = i+1; j < latticeNeighbors.size(); j++)
        {
            distances.push_back(roundDouble(structure.getDistance2(latticeNeighbors[i].index,
                                                                   latticeNeighbors[i].unitcellOffset,
                                                                   latticeNeighbors[j].index,
                                                                   latticeNeighbors[j].unitcellOffset)));
        }
    }
    Cluster cluster = Cluster(sites, distances);
    //cluster.print();
    _clusterCounts[cluster][elements] += 1;
}