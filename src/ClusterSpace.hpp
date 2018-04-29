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
This is the cluster space object.

It will have the definition of the cluster space a cluster expansion is based on.

*/

class ClusterSpace
{
  public:

    /// Constructor.
    ClusterSpace(std::vector<int>, std::vector<std::string>, const OrbitList);

    /// Returns the cluster vector corresponding to the input structure.
    std::vector<double> getClusterVector(const Structure &) const;

    /// Returns the full cluster product of entire cluster (elements vector). Assuming all sites have the same number of allowed components.
    double getClusterProduct(const std::vector<int> &mcVector,
                             const std::vector<int> &numberOfAllowedComponents,
                             const std::vector<int> &elements) const;

    ///setup  _clusterSpaceInfo
    void setupClusterSpaceInfo();

    /// Returns cluster space information.
    std::pair<int, std::vector<int>> getClusterSpaceInfo(const unsigned int);

    /// Returns the cluster space size, i.e. the length of a cluster vector.
    size_t getClusterSpaceSize();

    /// Returns the cutoff for each order.
    std::vector<double> getCutoffs() const { return _clusterCutoffs; }

    /// Returns list of elements associated with cluster space as chemical symbols.
    std::vector<std::string> getChemicalSymbols() const
    {
        std::vector<std::string> elements;
        for (const auto &intEl : _elements)
        {
            elements.push_back(PeriodicTable::intStr[intEl]);
        }
        return elements;
    }

    ///returns a orbit from the orbit list
    Orbit getOrbit(const size_t index) const
    {
        return _primitiveOrbitList.getOrbit(index);
    }

    OrbitList getOrbitList() const
    {
        return _primitiveOrbitList;
    }

    ///returns the primitive structure
    Structure getPrimitiveStructure() const
    {
        return _primitiveStructure;
    }

    ClusterCounts getNativeClusters(const Structure &structure) const;

    /// Return the MC vector permutations for each mc vector in the set of mc vectors.
    std::vector<std::vector<std::vector<int>>> getMCVectorPermutations(const std::vector<std::vector<int>> &, const int ) const;

    ///Returns the allowed occupations on the sites
    std::vector<int> getAllowedOccupations(const Structure &structure, const std::vector<LatticeSite> &latticeNeighbors) const;


  private:

    /// Sets up a map between chemical elements and the internal species enumeration scheme.
    void setupElementMap(const std::vector<std::string>);

    /// Elements identified by atomic number considered in this cluster space.
    std::vector<int> _elements;

    /// Maps atomic numbers to the internal species enumeration scheme.
    std::map<int, int> _elementMap;


  private:
    /// Number of allowed components on each site of the primitive structure.
    std::vector<int> _numberOfAllowedComponents;

    /// Primitive structure.
    Structure _primitiveStructure;

    /// Orbit list for primitive structure.
    OrbitList _primitiveOrbitList;

    /// Radial cutoffs by cluster order starting with pairs.
    std::vector<double> _clusterCutoffs;

    /// Cluster space information. The first index (int) corresponds to the orbit index, the second index (vector of ints) corresponds to multi-component vectors.
    std::vector<std::pair<int, std::vector<int>>> _clusterSpaceInfo;

    ///a boolean to keep track if this is initialized or not
    bool _isClusterSpaceInitialized;

    /**
    This maps a orbit to other orbits i.e. of _equalSitesList[0] = [1,2,3] then orbit zero should be seen as equal to orbit 1,2 and three.
    This will mean that when retrieving the cluster vector the zeroth element will be a combination of orbit 0,1,2 and 3.
    */
    std::map<int, std::vector<int>> _equalOrbitsList;

    ///return the default clusterfunction of the input element when this site can have Mi elements
    double defaultClusterFunction(const int Mi, const int clusterFunction, const int element) const;

};
