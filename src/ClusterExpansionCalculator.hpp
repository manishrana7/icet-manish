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
#include "ClusterCounts.hpp"
#include "PeriodicTable.hpp"
#include "Vector3dCompare.hpp"
using namespace Eigen;


class ClusterExpansionCalculator
{
    public:
    ClusterExpansionCalculator(const ClusterSpace &, const Structure &);

    // std::vector<double> getLocalClusterVector(const Structure &, const int);
    std::vector<double> getLocalClusterVector(const std::vector<int>& ,const int, const std::vector<int>, bool);
    double getLocalContribution(const Structure &, const int) const;
    OrbitList getLocalOrbitList(int index);
    void testRemovingSites();
    void checkNoSelfInteractions();

    private:
    std::map<Vector3d, OrbitList, Vector3dCompare> _localOrbitlists;

    void validateBasisAtomOrbitLists();
    ClusterSpace _clusterSpace;
    Structure _superCell;
    LocalOrbitListGenerator _theLog;
    OrbitList _fullPrimitiveOrbitList;
    ///this maps a latticeNeighbor from the primitive and get the equivalent in supercell
    std::unordered_map<LatticeSite, LatticeSite> _primToSupercellMap;

};