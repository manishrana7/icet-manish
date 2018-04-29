#include "Geometry.hpp"

namespace icet {

    /**
    @details The radius is defined as the maximum average distance to the center of mass of the cluster.
    @param latticeSites a list of lattice sites
    @param structure atomic configuration, which allows transforming LatticeSite information to a Cartesian coordinate
    */
    double getClusterRadius(const std::vector<LatticeSite> &latticeSites, const Structure &structure)
    {
        // compute the center of the cluster of lattice sites
        Vector3d centerPosition = {0.0, 0.0, 0.0};
        for(const auto &latnbr : latticeSites)
        {
            centerPosition += structure.getPosition(latnbr) / latticeSites.size();
        }

        // compute the average distance to the center
        double avgDistanceToCenter = 0.0;
        for(const auto &latnbr : latticeSites)
        {
            avgDistanceToCenter += (centerPosition - structure.getPosition(latnbr)).norm() / latticeSites.size();
        }
        return avgDistanceToCenter;
    }

}
