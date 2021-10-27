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
    Orbit(const std::vector<std::vector<LatticeSite>> &, const Structure &, const std::set<std::vector<int>> &);

    /// Adds one cluster to the orbit.
    void addEquivalentCluster(const std::vector<LatticeSite> &, bool = false);

    /// Returns the number of equivalent clusters in this orbit.
    size_t size() const { return _equivalentClusters.size(); }

    /// Returns the radius of the representative cluster in this orbit.
    double radius() const { return _representativeCluster.radius(); }

    /// Returns the sorted, representative cluster for this orbit.
    const Cluster &getRepresentativeCluster() const { return _representativeCluster; }

    /// Returns the sites that define the representative cluster of this orbit.
    const std::vector<LatticeSite> &getSitesOfRepresentativeCluster() const { return _equivalentClusters[0]; }

    /// Returns the equivalent cluster.
    std::vector<LatticeSite> getClusterByIndex(unsigned int) const;

    /// Returns all equivalent clusters.
    const std::vector<std::vector<LatticeSite>> &getEquivalentClusters() const { return _equivalentClusters; }

    /// Sets the equivalent clusters.
    void setEquivalentClusters(const std::vector<std::vector<LatticeSite>> &equivalentClusters) { _equivalentClusters = equivalentClusters; }

    /// Sorts equivalent clusters.
    void sort() { std::sort(_equivalentClusters.begin(), _equivalentClusters.end()); }

    /// Returns the number of bodies of the cluster that represent this orbit.
    unsigned int order() const { return _representativeCluster.order(); }

    /// Gets the allowed permutations of clusters.
    const std::set<std::vector<int>> &getAllowedClusterPermutations() const { return _allowedClusterPermutations; }

    /// Returns the relevant multicomponent vectors of this orbit given the number of allowed components.
    std::vector<std::vector<int>> getMultiComponentVectors(const std::vector<int> &Mi_local) const;

    std::vector<std::vector<int>> getAllPossibleMultiComponentVectorPermutations(const std::vector<int> &Mi_local) const;

    /// Returns true if the input sites exists in _equivalentClusters, order does not matter if sorted=false.
    bool contains(const std::vector<LatticeSite>, bool) const;

    /// Counts occupations of clusters in this orbit
    std::map<std::vector<int>, double> countClusters(const Structure &, int doNotDoubleCountThisSiteIndex = -1) const;

    /// Counts changes in the occupation of clusters in this orbit
    std::map<std::vector<int>, double> countClusterChanges(const Structure &, const int, const int) const;

    /// Comparison operator for automatic sorting in containers.
    friend bool operator==(const Orbit &orbit1, const Orbit &orbit2)
    {
        throw std::logic_error("Reached equal operator in Orbit");
    }

    /// Comparison operator for automatic sorting in containers.
    friend bool operator<(const Orbit &orbit1, const Orbit &orbit2)
    {
        throw std::logic_error("Reached < operator in Orbit");
    }

    /**
    Creates a copy of this orbit and translates all LatticeSite offsets in equivalent sites.
    This will also transfer any existing permutations directly, which should be fine since an offset does not change permutations to the prototype sites.
    */
    friend Orbit operator+(const Orbit &orbit, const Eigen::Vector3d &offset)
    {
        Orbit orbitOffset = orbit;
        for (auto &latticeSites : orbitOffset._equivalentClusters)
        {
            for (auto &latticeSite : latticeSites)
            {
                latticeSite = latticeSite + offset;
            }
        }
        return orbitOffset;
    }

    /// Appends an orbit to this orbit.
    Orbit &operator+=(const Orbit &orbit_rhs)
    {
        // Get representative sites
        auto rep_sites_rhs = orbit_rhs.getSitesOfRepresentativeCluster();
        auto rep_sites_this = getSitesOfRepresentativeCluster();

        if (rep_sites_this.size() != rep_sites_rhs.size())
        {
            throw std::runtime_error("Orbit order is not equal (Orbit &operator+=)");
        }

        const auto rhsEquivalentClusters = orbit_rhs.getEquivalentClusters();

        // Insert rhs eq sites and corresponding permutations
        _equivalentClusters.insert(_equivalentClusters.end(), rhsEquivalentClusters.begin(), rhsEquivalentClusters.end());
        return *this;
    }

private:
    /// Container of equivalent sites for this orbit
    std::vector<std::vector<LatticeSite>> _equivalentClusters;

    /// Representative sorted cluster for this orbit
    Cluster _representativeCluster;

    /// Contains the allowed sites permutations. i.e. if 0,2,1 is in this set then 0,1,0 is the same MC vector as 0,0,1
    std::set<std::vector<int>> _allowedClusterPermutations;

    /// Check whether a site is included in a vector of lattice sites
    bool isSiteIncluded(const int, const std::vector<LatticeSite> &) const;
};

namespace std
{
    /// Stream operator.
    ostream &operator<<(ostream &, const Orbit &);
}
