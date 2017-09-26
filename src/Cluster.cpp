#include "Cluster.hpp"


Cluster::Cluster(const Structure &structure, const std::vector<LatticeNeighbor> &latticeNeighbors,
        const bool sortedCluster, const int clusterTag )
{
    _symprec = 1e-6;
    size_t clusterSize = latticeNeighbors.size();
    std::vector<int> sites(clusterSize);
    std::vector<double> distances;
    distances.reserve((clusterSize * (clusterSize - 1) / 2));
    Vector3d avg_position = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < latticeNeighbors.size(); i++)
    {
        avg_position += structure.getPosition(latticeNeighbors[i]);
        sites[i] = structure.getSite(latticeNeighbors[i].index);
        for (size_t j = i + 1; j < latticeNeighbors.size(); j++)
        {
            double distance = roundDouble(structure.getDistance2(latticeNeighbors[i].index,
                                                                   latticeNeighbors[i].unitcellOffset,
                                                                   latticeNeighbors[j].index,
                                                                   latticeNeighbors[j].unitcellOffset));
            
            distances.push_back(distance);
        }
    }
    //from avg_position (center position ) get average distance to center
    double meanDistanceToCenter = 0.0;
    for(const auto &latNbr :latticeNeighbors)
    {
        meanDistanceToCenter += (structure.getPosition(latNbr) - avg_position).norm();
    }

    meanDistanceToCenter /= (double) latticeNeighbors.size();

    _sites = sites;
    _distances = distances;
    _sortedCluster = sortedCluster;
    _clusterTag = clusterTag;
    _geometricalSize = meanDistanceToCenter;
    if (_sortedCluster)
    {
        sortCluster();
    }
}