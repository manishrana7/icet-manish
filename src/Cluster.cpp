#include "Cluster.hpp"

/**
@details Create an instance of a cluster.
@param structure icet structure object
@param latticeSites list of lattice sites that form the cluster
*/
Cluster::Cluster(const Structure &structure,
                 const std::vector<LatticeSite> &latticeSites)
{
    _latticeSites = latticeSites;
    //_structure = structure;

    std::vector<double> distances;
    float thisOrder = order();
    distances.reserve((thisOrder * (thisOrder - 1) / 2));
    for (size_t i = 0; i < thisOrder; i++)
    {
        for (size_t j = i + 1; j < thisOrder; j++)
        {
            double distance = structure.getDistance(_latticeSites[i].index(),
                                                    _latticeSites[j].index(),
                                                    _latticeSites[i].unitcellOffset(),
                                                    _latticeSites[j].unitcellOffset());
            distances.push_back(distance);
        }
    }
    _distances = distances;

    _radius = icet::getGeometricalRadius(_latticeSites, structure);
}

std::vector<double> Cluster::distances() const
{
    return _distances;
    /*
    std::vector<double> distances;
    float thisOrder = order();
    distances.reserve((thisOrder * (thisOrder - 1) / 2));
    for (size_t i = 0; i < thisOrder; i++)
    {
        for (size_t j = i + 1; j < thisOrder; j++)
        {
            double distance = _structure.getDistance(_latticeSites[i].index(),
                                                     _latticeSites[j].index(),
                                                     _latticeSites[i].unitcellOffset(),
                                                     _latticeSites[j].unitcellOffset());
            distances.push_back(distance);
        }
    }
    return distances;
    */
}

double Cluster::radius() const
{
    return _radius;
    //return icet::getGeometricalRadius(_latticeSites, _structure);
}

namespace std
{
    /// Stream operator.
    ostream &operator<<(ostream &os, const Cluster &cluster)
    {
        for (const auto d : cluster.distances())
        {
            os << d << " ";
        }
        os << cluster.radius();
        return os;
    }
}
