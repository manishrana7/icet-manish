#include "ClusterCounts.hpp"

///
void ClusterCounts::count_using_mbnl(const Structure &structure, ManybodyNeighborlist &mbnl, const int order)
{
    bool bothways = false;

    for (size_t latticeIndex = 0; latticeIndex < structure.size(); latticeIndex++)
    {
        auto latticeNeighbors = mbnl.build(structure, latticeIndex, order, bothways);
        void countLatticeNeighbors(structure, latticeNeighbors);
    }
}

/**

Counts clusters given this compact form of latticeneighbors (see ManybodyNeighborlist for more details)
*/ 
// build(const Neighborlist &nl, int index, int order, bool);
void ClusterCounts::countLatticeNeighbors(const Structure &structure,
                                          const std::vector<std::pair<std::vector<std::pair<int, Vector3d>>, std::vector<std::pair<int, Vector3d>>>> &latticeNeighbors)
{
    for (const auto &neighborPair : latticeNeighbors)
    {
        //Now we have std::pair<std::vector<std::pair<int, Vector3d>>, std::vector<std::pair<int, Vector3d>>>
        //pair.first == the base indices and pair.second is all indices that form clusters with the base indices
        for (const auto &combinationIndice : neighborPair.second)
        {
            auto latticePointsForCluster = neighborPair.second;
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
                          const std::vector<std::pair<int, Vector3d>> &latticeNeighbors)
{

Cluster cluster;

}