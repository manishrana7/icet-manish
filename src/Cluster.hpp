#pragma once

#include <boost/functional/hash.hpp>
#include "FloatType.hpp"
#include "Geometry.hpp"
#include "LatticeSite.hpp"
#include "Structure.hpp"

using boost::hash;
using boost::hash_combine;

/// This class handles information pertaining to a single cluster.
class Cluster
{
public:
    /// Empty constructor.
    Cluster() {}

    /// Creates cluster from a structure and a set of lattice sites.
    Cluster(const std::shared_ptr<Structure>,
            const std::vector<LatticeSite> &);

    /// Returns the lattice sites of this cluster.
    std::vector<LatticeSite> getLatticeSites() const { return _latticeSites; }

    /// Returns the order (i.e., the number of sites) of the cluster.
    unsigned int order() const { return _latticeSites.size(); }

    /// Returns the radius of the cluster.
    double radius() const;

    /// Returns distances between points in the cluster.
    std::vector<double> distances() const;

    /// Comparison operator for automatic sorting.
    friend bool operator<(const Cluster &cluster1, const Cluster &cluster2)
    {
        return cluster1.getLatticeSites() < cluster2.getLatticeSites();
    }

private:
    /// The lattice sites in the cluster.
    std::vector<LatticeSite> _latticeSites;

    /// The structure that the lattice sites refer to.
    std::shared_ptr<Structure> _structurePtr;

    std::vector<double> _distances;

    float _radius;
};
