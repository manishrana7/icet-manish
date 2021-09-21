#pragma once

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <Eigen/Dense>

#include "Cluster.hpp"
#include "LatticeSite.hpp"
//#include "ManyBodyNeighborList.hpp"
#include "OrbitList.hpp"
//#include "PeriodicTable.hpp"
#include "Structure.hpp"

using namespace Eigen;

namespace py = pybind11;

/// This class is used to generate a count of the number of different clusters.
class ClusterCounts
{
public:
  /// Constructor.
  ClusterCounts() {}
  void count(const Structure &, const std::vector<std::vector<LatticeSite>> &, const Cluster &, bool, int doNotDoubleCountThisSiteIndex = -1);
  void countChange(const Structure &, const int, const int, const std::vector<std::vector<LatticeSite>> &, const Cluster &, bool, int doNotDoubleCountThisSiteIndex = -1);
  void countOrbitList(const Structure &, const OrbitList &, bool keepOrder, bool permuteSites = false, int maxOrbit = -1, int doNotDoubleCountThisSiteIndex = -1);
  void countOrbitListChange(const Structure &, const int, const int, const OrbitList &, bool keepOrder, int maxOrbit = -1, int doNotDoubleCountThisSiteIndex = -1);

  /**
    @details Returns a map representing the cluster counts. The key of the map
              represents a cluster, while the value is another map. In the
              latter the key represents the chemical species, while the value
              is the cluster count.
    **/
  std::unordered_map<Cluster, std::map<std::vector<int>, double>> getClusterCounts() const { return _clusterCounts; }

  /// Returns the cluster counts size i.e. the total number of clusters
  size_t size() const
  {
    return _clusterCounts.size();
  }

  /// Reset cluster counts
  void reset()
  {
    _clusterCounts.clear();
  }

  /**
    @details Map representing cluster counts. The key represents a cluster,
              while the value is another map. In the latter the key represents
              the chemical species, while the value is the cluster count.
    @todo make private
    **/
  std::unordered_map<Cluster, std::map<std::vector<int>, double>> _clusterCounts;
};
