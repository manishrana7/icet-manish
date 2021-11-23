#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "Cluster.hpp"
#include "LatticeSite.hpp"
#include "Symmetry.hpp"
#include "VectorOperations.hpp"

using namespace Eigen;

/**
This class handles an orbit.

An orbit is a set of clusters that are equivalent under the symmetry operations
of the underlying lattice. Each cluster is represented by a set of lattice
sites. An orbit is characterized by a representative cluster.

*/

class Orbit
{
public:
    /// Constructor.
    Orbit(const std::vector<Cluster>, const std::set<std::vector<int>>);

    /// Adds one cluster to the orbit.
    void addCluster(const Cluster &);

    /// Returns the number of clusters in this orbit.
    size_t size() const { return _clusters.size(); }

    /// Returns the radius of the representative cluster in this orbit.
    double radius() const { return getRepresentativeCluster().radius(); }

    /// Returns the representative cluster for this orbit
    const Cluster &getRepresentativeCluster() const { return _representativeCluster; }

    /// Returns all clusters in this orbit.
    const std::vector<Cluster> &getClusters() const { return _clusters; }

    /// Returns the number of bodies of the cluster that represent this orbit.
    unsigned int order() const { return getRepresentativeCluster().order(); }

    /// Gets the allowed permutations of clusters.
    std::set<std::vector<int>> getAllowedClusterPermutations() const { return _allowedClusterPermutations; }

    /// Returns the relevant multicomponent vectors of this orbit given the number of allowed components.
    std::vector<std::vector<int>> getMultiComponentVectors(const std::vector<int> &Mi_local) const;

    std::vector<std::vector<int>> getAllPossibleMultiComponentVectorPermutations(const std::vector<int> &Mi_local) const;

    /// Returns true if the input sites exists in _clusters, order does not matter if sorted=false.
    bool contains(const std::vector<LatticeSite>, bool) const;

    /// Counts occupations of clusters in this orbit.
    std::map<std::vector<int>, double> getClusterCounts(std::shared_ptr<Structure>, int doNotDoubleCountThisSiteIndex = -1) const;

    /// Counts changes in the occupation of clusters in this orbit.
    std::map<std::vector<int>, double> getClusterCountChanges(std::shared_ptr<Structure>, const int, const int) const;

    /// Returns a copy of this orbit in the given (supercell) structure.
    void transformToSupercell(std::shared_ptr<Structure>,
                              std::unordered_map<LatticeSite, LatticeSite> &,
                              const double);

    /// Translates the orbit with an offset
    void translate(const Vector3d &);

    /// Comparison operator for automatic sorting in containers.
    friend bool
    operator==(const Orbit &orbit1, const Orbit &orbit2)
    {
        throw std::logic_error("Reached equal operator in Orbit");
    }

    /// Comparison operator for automatic sorting in containers.
    friend bool operator<(const Orbit &orbit1, const Orbit &orbit2)
    {
        throw std::logic_error("Reached < operator in Orbit");
    }

    /// Appends an orbit to this orbit.
    Orbit &operator+=(const Orbit &orbit_rhs)
    {
        // Get representative sites
        auto rep_sites_rhs = orbit_rhs.getRepresentativeCluster().getLatticeSites();
        auto rep_sites_this = _representativeCluster.getLatticeSites();

        if (rep_sites_this.size() != rep_sites_rhs.size())
        {
            throw std::runtime_error("Orbit order is not equal (Orbit &operator+=)");
        }

        const auto rhsClusters = orbit_rhs.getClusters();

        // Insert rhs eq sites
        _clusters.insert(_clusters.end(), rhsClusters.begin(), rhsClusters.end());
        return *this;
    }

private:
    /// Container of all clusters in this orbit
    std::vector<Cluster> _clusters;

    /// One of the clusters chosen to represent the orbit
    Cluster _representativeCluster;

    /// Contains the allowed sites permutations. i.e. if 0, 2, 1 is in this set then 0, 1, 0 is the same multi-component vector as 0, 0, 1
    std::set<std::vector<int>> _allowedClusterPermutations;
};

namespace std
{
    /// Stream operator.
    ostream &operator<<(ostream &, const Orbit &);
}
