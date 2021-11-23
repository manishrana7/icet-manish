#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include "Structure.hpp"
#include "ClusterSpace.hpp"
#include "OrbitList.hpp"
#include "LocalOrbitListGenerator.hpp"
#include "PeriodicTable.hpp"
#include "VectorOperations.hpp"
using namespace Eigen;

/**
@details This class provides a cluster expansion calculator. A cluster
    expansion calculator is specific for a certain supercell. Upon
    initialization various quantities specific to the given supercell are
    precomputed. This greatly speeds up subsequent calculations and enables one
    to carry out e.g., Monte Carlo simulations in a computationally efficient
    manner.
**/
class ClusterExpansionCalculator
{
public:
    /// Constructor.
    ClusterExpansionCalculator(const ClusterSpace &, const Structure &, const double);

    /// Returns change in cluster vector upon flipping occupation of one site
    std::vector<double> getClusterVectorChange(const py::array_t<int> &, const int, const int);

    /// Returns the full cluster vector.
    std::vector<double> getClusterVector(const py::array_t<int> &);

    /// Returns a local cluster vector; the contribution to the cluster vector from one site.
    std::vector<double> getLocalClusterVector(const py::array_t<int> &, int);

private:
    /// Maps offsets to local orbit lists.
    std::unordered_map<Vector3d, OrbitList, Vector3dHash> _localOrbitlists;

    /// Internal cluster space.
    ClusterSpace _clusterSpace;

    /// The supercell the calculator is created for.
    std::shared_ptr<Structure> _supercell;

    /// The full primitive orbit list, contains all clusters for the primitive cell.
    OrbitList _fullPrimitiveOrbitList;

    /// Maps a lattice site from the primitive and get the equivalent in the supercell.
    std::unordered_map<LatticeSite, LatticeSite> _primitiveToSupercellMap;

    /// Maps supercell index to its corresponding primitive cell offset.
    std::map<int, Vector3d> _indexToOffset;

    /// Placeholder for translated orbitlist
    OrbitList _translatedOrbitList;

    /// The full orbit list used when calculating full cluster vector
    OrbitList _fullOrbitList;
};
