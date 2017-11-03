#include "Geometry.hpp"


namespace icet
{
/// Returns the geomtetrical radius from the vectors of latticeneighbors and a unitcell
double getGeometricalRadius(const std::vector<LatticeNeighbor> &latticeNeigbhors, const Structure &structure)
{
    Vector3d centerPosition = {0.0, 0.0, 0.0};
    for(const auto &latnbr : latticeNeigbhors)
    {
        centerPosition += structure.getPosition(latnbr)/latticeNeigbhors.size();
    }

    double avgDistanceToCenter = 0.0;
    for(const auto &latnbr : latticeNeigbhors)
    {
        avgDistanceToCenter += (centerPosition-structure.getPosition(latnbr)).norm()/latticeNeigbhors.size();
    }
    return avgDistanceToCenter;
}



}
