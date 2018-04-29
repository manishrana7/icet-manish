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
#include "OrbitList.hpp"
#include "LocalOrbitListGenerator.hpp"
#include "ClusterCounts.hpp"
#include "PeriodicTable.hpp"

using namespace Eigen;

/**
@brief This class handles the cluster space.
@details It provides functionality for setting up a cluster space, calculating cluster vectors as well as retrieving various types of associated information.
*/

class ClusterSpace
{
  public:

    /// Constructor.
    ClusterSpace(std::vector<int>, std::vector<std::string>, const OrbitList);

    /// Returns the cluster vector corresponding to the input structure.
    std::vector<double> getClusterVector(const Structure &) const;

    /// Returns information concerning the cluster space.
    std::pair<int, std::vector<int>> getClusterSpaceInfo(const unsigned int);

    /// Returns the entire orbit list.
    OrbitList getOrbitList() const { return _primitiveOrbitList; }

    /// Returns an orbit from the orbit list.
    Orbit getOrbit(const size_t index) const { return _primitiveOrbitList.getOrbit(index); }

    /// Returns the native clusters.
    /// @todo What is a native cluster?
    ClusterCounts getNativeClusters(const Structure &structure) const;

    /// Returns the multi-component (MC) vector permutations for each MC vector in the set of MC vectors.
    /// @todo Clean up this description.
    std::vector<std::vector<std::vector<int>>> getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &, const int ) const;


  public:

    /// Returns the cutoff for each order.
    std::vector<double> getCutoffs() const { return _clusterCutoffs; }

    /// Returns the primitive structure.
    Structure getPrimitiveStructure() const { return _primitiveStructure; }

    /// Returns the number of allowed components for each site.
    std::vector<int> getNumberOfAllowedComponentsBySite(const Structure &structure, const std::vector<LatticeSite> &latticeNeighbors) const;

    /// Returns a list of elements associated with cluster space as chemical symbols.
    std::vector<std::string> getChemicalSymbols() const
    {
        std::vector<std::string> elements;
        for (const auto &intEl : _elements)
        {
            elements.push_back(PeriodicTable::intStr[intEl]);
        }
        return elements;
    }

    /// Returns the cluster space size, i.e. the length of a cluster vector.
    size_t getClusterSpaceSize()
    {
        if (!_isClusterSpaceInitialized)
        {
            collectClusterSpaceInfo();
        }
        return _clusterSpaceInfo.size();
    }

    /// Returns the mapping between atomic numbers and the internal species enumeration scheme.
    std::map<int, int> getElementMap() const
    {
        return _elementMap;
    }


  private:

    /// Returns the cluster product.
    /// @todo What is a cluster product?
    double getClusterProduct(const std::vector<int> &mcVector, const std::vector<int> &numberOfAllowedComponents, const std::vector<int> &elements) const;

    /// Collect information about the cluster space.
    void collectClusterSpaceInfo();

    /// Sets up a map between chemical elements and the internal species enumeration scheme.
    void setupElementMap(const std::vector<std::string>);

    /// Returns the default cluster function.
    double defaultClusterFunction(const int numberOfAllowedComponents, const int clusterFunction, const int element) const;


  private:

    /// True if cluster space has been initialized.
    bool _isClusterSpaceInitialized = false;

    /// Cluster space information.
    /// The first index (int) corresponds to the orbit index, the second index (vector of ints) refers to a multi-component vector.
    /// @todo Check description.
    std::vector<std::pair<int, std::vector<int>>> _clusterSpaceInfo;

    /// Number of allowed components on each site of the primitive structure.
    std::vector<int> _numberOfAllowedComponents;

    /// Primitive structure.
    Structure _primitiveStructure;

    /// List of orbits for primitive structure.
    OrbitList _primitiveOrbitList;

    /// Radial cutoffs by cluster order starting with pairs.
    std::vector<double> _clusterCutoffs;

    /// Elements identified by atomic number considered in this cluster space.
    std::vector<int> _elements;

    /// Map between atomic numbers and the internal species enumeration scheme.
    std::map<int, int> _elementMap;

};
