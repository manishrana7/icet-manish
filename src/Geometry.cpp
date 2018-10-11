#include "Geometry.hpp"

namespace icet
{

    /**
     @details This function computes the geometrical radius of a cluster that is defined in terms of lattice sites.
    @param latticeSites a list of lattice sites
    @param structure atomic configuration used to convert information in latticeSites to Cartesian coordinates
    */
    double getGeometricalRadius(const std::vector<LatticeSite> &latticeSites, const Structure &structure)
    {
        // Compute the center of the cluster.
        Vector3d centerPosition = {0.0, 0.0, 0.0};
        for(const auto &latnbr : latticeSites)
        {
            centerPosition += structure.getPosition(latnbr) / latticeSites.size();
        }

        // Compute the average distance of the points in the cluster to its center.
        double avgDistanceToCenter = 0.0;
        for(const auto &latnbr : latticeSites)
        {
            avgDistanceToCenter += (centerPosition - structure.getPosition(latnbr)).norm() / latticeSites.size();
        }
        return avgDistanceToCenter;
    }

}
