#pragma once

#include "Structure.hpp"
#include "OrbitList.hpp"
#include "LocalOrbitListGenerator.hpp"
#include "ClusterCounts.hpp"
#include "PeriodicTable.hpp"

//namespace icet {

using namespace std;

/**
@brief This class handles the cluster space.
@details It provides functionality for setting up a cluster space, calculating cluster vectors as well as retrieving various types of associated information.
@todo Add mathematical definitions.
*/

class ClusterSpace
{
  public:

    /// Constructor.
    ClusterSpace(vector<int>, vector<string>, const OrbitList);

    /// Returns the cluster vector corresponding to the input structure.
    vector<double> getClusterVector(const Structure &) const;

    /// Returns information concerning the cluster space.
    pair<int, vector<int>> getClusterSpaceInfo(const unsigned int);

    /// Returns the entire orbit list.
    OrbitList getOrbitList() const { return _orbitList; }

    /// Returns an orbit from the orbit list.
    Orbit getOrbit(const size_t index) const { return _orbitList.getOrbit(index); }

    /// Returns the native clusters.
    /// @todo What is a native cluster? Partial answer: clusters within the unit cell?
    ClusterCounts getNativeClusters(const Structure &structure) const;

    /// Returns the multi-component (MC) vector permutations for each MC vector in the set of MC vectors.
    /// @todo Clean up this description.
    vector<vector<vector<int>>> getMultiComponentVectorPermutations(const vector<vector<int>> &, const int ) const;


  public:

    /// Returns the cutoff for each order.
    vector<double> getCutoffs() const { return _clusterCutoffs; }

    /// Returns the primitive structure.
    Structure getPrimitiveStructure() const { return _primitiveStructure; }

    /// Returns the number of allowed components for each site.
    vector<int> getNumberOfAllowedComponentsForEachSite(const Structure &structure, const vector<LatticeSite> &latticeNeighbors) const;

    /// Returns a list of elements associated with cluster space as chemical symbols.
    vector<string> getChemicalSymbols() const
    {
        vector<string> elements;
        for (const auto &elem : _elements)
            elements.push_back(PeriodicTable::intStr[elem]);
        return elements;
    }

    /// Returns the cluster space size, i.e. the length of a cluster vector.
    size_t getClusterSpaceSize()
    {
        if (!_isClusterSpaceInitialized)
            collectClusterSpaceInfo();
        return _clusterSpaceInfo.size();
    }

    /// Returns the mapping between atomic numbers and the internal species enumeration scheme.
    map<int, int> getElementMap() const { return _elementMap; }


  private:

    /// Returns the cluster product.
    /// @todo What is a cluster product?
    double getClusterProduct(const vector<int> &mcVector, const vector<int> &numberOfAllowedComponents, const vector<int> &elements) const;

    /// Collect information about the cluster space.
    void collectClusterSpaceInfo();

    /// Sets up a map between chemical elements and the internal species enumeration scheme.
    void setupElementMap(const vector<string>);

    /// Returns the default cluster function.
    double defaultClusterFunction(const int numberOfAllowedComponents, const int clusterFunction, const int element) const;


  private:

    /// True if cluster space has been initialized.
    bool _isClusterSpaceInitialized = false;

    /// Cluster space information.
    /// The first index (int) corresponds to the orbit index, the second index (vector of ints) refers to a multi-component vector.
    /// @todo Check description.
    /// @todo This function returns a very specific type of information. Consider giving it a more descriptive name.
    vector<pair<int, vector<int>>> _clusterSpaceInfo;

    /// Number of allowed components on each site of the primitive structure.
    vector<int> _numberOfAllowedComponents;

    /// Primitive (prototype) structure.
    Structure _primitiveStructure;

    /// List of orbits based on primitive structure.
    OrbitList _orbitList;

    /// Radial cutoffs by cluster order starting with pairs.
    vector<double> _clusterCutoffs;

    /// Elements identified by atomic number considered in this cluster space.
    vector<int> _elements;

    /// Map between atomic numbers and the internal species enumeration scheme.
    map<int, int> _elementMap;

};

//}