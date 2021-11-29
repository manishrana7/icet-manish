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
@brief This struct keeps track of information pertaining to a specific element
       in the cluster vector.
*/
struct ClusterVectorElement
{
    /// A multi-component vector contains the indices of the point functions
    /// (only non-trivial if the number of components are more than two)
    std::vector<int> multicomponentVector;

    /// Site permutations describe how the sites in the cluster can be re-ordered
    std::vector<std::vector<int>> sitePermutations;

    /// Multiplicity for this cluster vector element
    size_t multiplicity;
};

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
    double radius() const { return representativeCluster().radius(); }

    /// Returns the representative cluster for this orbit
    const Cluster &representativeCluster() const { return _representativeCluster; }

    /// Returns all clusters in this orbit.
    const std::vector<Cluster> &clusters() const { return _clusters; }

    /// Returns the number of bodies of the cluster that represent this orbit.
    unsigned int order() const { return representativeCluster().order(); }

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
    void translate(const Vector3i &);

    /// Appends an orbit to this orbit.
    Orbit &operator+=(const Orbit &orbit_rhs);

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

    const std::vector<ClusterVectorElement> &clusterVectorElements() const { return _clusterVectorElements; }

    std::vector<ClusterVectorElement> _clusterVectorElements;

private:
    /// Container of all clusters in this orbit
    std::vector<Cluster> _clusters;

    /// One of the clusters chosen to represent the orbit
    Cluster _representativeCluster;

    /// Contains the allowed sites permutations. i.e. if 0, 2, 1 is in this set then 0, 1, 0 is the same multi-component vector as 0, 0, 1
    std::set<std::vector<int>> _allowedClusterPermutations;

    ///
    void _computeMultiComponentVectors();

    std::vector<std::vector<std::vector<int>>> _getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &) const;
};

namespace std
{
    /// Stream operator.
    ostream &operator<<(ostream &, const Orbit &);
}
