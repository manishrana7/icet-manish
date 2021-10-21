#include "Cluster.hpp"

/**
@details Create an instance of a cluster.
@param structure icet structure object
@param latticeSites list of lattice sites that form the cluster
*/
Cluster::Cluster(const Structure &structure,
                 const std::vector<LatticeSite> &latticeSites)
{
    _order = latticeSites.size();
    std::vector<double> distances;
    distances.reserve((_order * (_order - 1) / 2));
    for (size_t i = 0; i < latticeSites.size(); i++)
    {
        for (size_t j = i + 1; j < latticeSites.size(); j++)
        {
            double distance = structure.getDistance(latticeSites[i].index(),
                                                    latticeSites[j].index(),
                                                    latticeSites[i].unitcellOffset(),
                                                    latticeSites[j].unitcellOffset());
            distances.push_back(distance);
        }
    }
    _distances = distances;
    _radius = icet::getGeometricalRadius(latticeSites, structure);
}

namespace std
{
    /// Stream operator.
    ostream& operator<<(ostream& os, const Cluster& cluster)
    {
        for (const auto d : cluster.distances())
        {
            os << d << " ";
        }
        os << cluster.radius();
        return os;
    }
}
