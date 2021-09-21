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

/**
    @details Returns a map representing the cluster counts. The key of the map
              represents a cluster, while the value is another map. In the
              latter the key represents the chemical species, while the value
              is the cluster count.
    **/

/// This class is used to generate a count of the number of different clusters.
std::map<std::vector<int>, double> countClusters(const Structure &, const std::vector<std::vector<LatticeSite>> &, bool, int doNotDoubleCountThisSiteIndex = -1);
std::map<std::vector<int>, double> countClusterChanges(const Structure &, const int, const int, const std::vector<std::vector<LatticeSite>> &, bool, int doNotDoubleCountThisSiteIndex = -1);
std::unordered_map<Cluster, std::map<std::vector<int>, double>> countOrbitList(const Structure &, const OrbitList &, bool keepOrder, bool permuteSites = false, int maxOrbit = -1, int doNotDoubleCountThisSiteIndex = -1);
std::unordered_map<Cluster, std::map<std::vector<int>, double>> countOrbitListChange(const Structure &, const int, const int, const OrbitList &, bool keepOrder, int maxOrbit = -1, int doNotDoubleCountThisSiteIndex = -1);
