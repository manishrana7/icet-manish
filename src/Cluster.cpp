#include "Cluster.hpp"

/**
@details Create an instance of a cluster.
@param structure icet structure object
@param latticeSites list of lattice sites that form the cluster
*/
Cluster::Cluster(const std::vector<LatticeSite> &latticeSites,
                 std::shared_ptr<const Structure> structure)
{
    _latticeSites = latticeSites;
    _structure = structure;
}

/** 
@brief Translate this cluster by an offset.
@param offset Coordinates referring to the axes of the structure in this cluster
**/
void Cluster::translate(const Eigen::Vector3d &offset)
{
    for (LatticeSite &site : _latticeSites)
    {
        site.addUnitcellOffset(offset);
    }
}

/**
@details Transforms a site from the primitive structure to a given supercell.
This involves finding a map from the site in the primitive cell to the supercell.
If no map is found mapping is attempted based on the position of the site in the supercell.
(The map is important for performance.)
@param supercell supercell structure
@param primToSuperMap map from primitive to supercell
@param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates
**/
void Cluster::transformToSupercell(std::shared_ptr<const Structure> supercell,
                                   std::unordered_map<LatticeSite, LatticeSite> &primToSuperMap,
                                   const double fractionalPositionTolerance)
{
    LatticeSite supercellSite;
    Vector3d sitePosition;
    for (LatticeSite &site : _latticeSites)
    {
        auto find = primToSuperMap.find(site);

        if (find == primToSuperMap.end())
        {
            sitePosition = _structure->getPosition(site);
            supercellSite = supercell->findLatticeSiteByPosition(sitePosition, fractionalPositionTolerance);
            primToSuperMap[site] = supercellSite;
        }
        else
        {
            supercellSite = primToSuperMap[site];
        }

        // overwrite site to match supercell index offset
        site.setIndex(supercellSite.index());
        site.setUnitcellOffset(supercellSite.unitcellOffset());
    }
    // Now all sites refer to the supercell, so we change the _structure to point to the supercell
    _structure = supercell;
}

/**
@brief This function computes the geometrical radius of the cluster.
*/
double Cluster::radius() const
{
    // Compute the center of the cluster.
    Vector3d centerPosition = {0.0, 0.0, 0.0};
    for (const auto &site : _latticeSites)
    {   
        centerPosition += _structure->getPosition(site);
    }
    centerPosition /= _latticeSites.size();

    // Compute the average distance of the points in the cluster to its center.
    double avgDistanceToCenter = 0.0;
    for (const auto &site : _latticeSites)
    {
        avgDistanceToCenter += (centerPosition - _structure->getPosition(site)).norm();
    }
    avgDistanceToCenter /= _latticeSites.size();
    return avgDistanceToCenter;
}

/**
@brief Returns the positions of the sites in this cluster in Cartesian coordinates.
**/
std::vector<Vector3d> Cluster::getPositions() const {
    std::vector<Vector3d> positions;
    for (const LatticeSite &site : _latticeSites) {
        positions.push_back(_structure->getPosition(site));
    }
    return positions;
}

/**
@brief Check whether a site index is included with a zero offset.
@param siteIndex Index of site to check whether it is included
*/
bool Cluster::isSiteIndexIncludedWithZeroOffset(int siteIndex) const
{
    return std::any_of(_latticeSites.begin(), _latticeSites.end(), [=](const LatticeSite &ls)
                       { return ls.index() == siteIndex && ls.unitcellOffset().norm() < 1e-4; });
}

/**
@brief Count the number of occurences of a site index among the sites in this cluster
@param siteIndex Index of site to count
*/
unsigned int Cluster::countOccurencesOfSiteIndex(int siteIndex) const
{
    return std::count_if(_latticeSites.begin(), _latticeSites.end(), [=](const LatticeSite &ls)
                         { return ls.index() == siteIndex; });
}