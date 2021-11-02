#include "Cluster.hpp"

/**
@details Create an instance of a cluster.
@param structure icet structure object
@param latticeSites list of lattice sites that form the cluster
*/
Cluster::Cluster(const std::shared_ptr<Structure> structurePtr,
                 const std::vector<LatticeSite> &latticeSites)
{
    _latticeSites = latticeSites;
    _structurePtr = structurePtr;
}

std::vector<double> Cluster::distances() const
{
    std::vector<double> distances;
    float thisOrder = order();
    distances.reserve((thisOrder * (thisOrder - 1) / 2));
    for (size_t i = 0; i < thisOrder; i++)
    {
        for (size_t j = i + 1; j < thisOrder; j++)
        {
            double distance = (*_structurePtr).getDistance(_latticeSites[i].index(), _latticeSites[j].index(), _latticeSites[i].unitcellOffset(), _latticeSites[j].unitcellOffset());
            distances.push_back(distance);
        }
    }
    return distances;
}

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
@param supercell supercell structure
@param primToSuperMap map from primitive to supercell
@param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates
**/
void Cluster::transformSitesToSupercell(const Structure &supercell,
                                        std::unordered_map<LatticeSite, LatticeSite> &primToSuperMap,
                                        const double fractionalPositionTolerance)
{
    for (LatticeSite &site : _latticeSites)
    {
        auto find = primToSuperMap.find(site);
        LatticeSite supercellSite;
        if (find == primToSuperMap.end())
        {
            Vector3d sitePosition = (*_structurePtr).getPosition(site);
            supercellSite = supercell.findLatticeSiteByPosition(sitePosition, fractionalPositionTolerance);
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
}

double Cluster::radius() const
{
    return icet::getGeometricalRadius(_latticeSites, *_structurePtr);
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
