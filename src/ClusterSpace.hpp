#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "LocalOrbitListGenerator.hpp"
#include "OrbitList.hpp"
#include "PeriodicTable.hpp"
#include "Structure.hpp"
#include "VectorOperations.hpp"

/**
@brief
*/
struct ClusterVectorElementInfo
{
  std::vector<int> multiComponentVector;
  std::vector<std::vector<int>> sitePermutations;
  int clusterVectorIndex;
  double multiplicity;
};

/**
@brief This class handles the cluster space.
@details It provides functionality for setting up a cluster space, calculating
cluster vectors as well as retrieving various types of associated information.
*/

class ClusterSpace
{
public:
  /// Constructor.
  ClusterSpace(){};
  ClusterSpace(std::vector<std::vector<std::string>> &, const OrbitList &, const double, const double);

  /// Returns the cluster vector corresponding to the input structure.
  std::vector<double> getClusterVector(const Structure &, const double) const;

  /// Returns information concerning the association between orbits and multi-component vectors.
  std::pair<int, std::vector<int>> getMultiComponentVectorsByOrbit(const unsigned int);

  /// Returns the entire primitive orbit list.
  const OrbitList &getPrimitiveOrbitList() const { return _primitiveOrbitList; }

  /// Returns an orbit from the orbit list.
  const Orbit getOrbit(const size_t index) const { return _primitiveOrbitList.getOrbit(index); }

  /// Returns the native clusters.
  /// @todo What is a native cluster? Partial answer: clusters within the unit cell?
  std::unordered_map<Cluster, std::map<std::vector<int>, double>> getNativeClusters(const Structure &structure) const;

  /// Returns the multi-component (MC) vector permutations for each MC vector in the set of input vectors.
  /// @todo Clean up this description.
  std::vector<std::vector<std::vector<int>>> getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &, const int) const;

  /// Returns the cluster vector given the orbit list and a structure
  const std::vector<double> occupyClusterVector(const OrbitList &, const Structure &, const double firstElement = 1.0, const int flipIndex = -1, const int newOccupation = -1, const bool permuteClusters = false) const;

  /// Returns the cutoff for each order.
  std::vector<double> getCutoffs() const { return _clusterCutoffs; }

  /// Returns the primitive structure.
  const Structure &getPrimitiveStructure() const { return _primitiveStructure; }

  /// Returns the number of allowed components for each site.
  std::vector<int> getNumberOfAllowedSpeciesBySite(const Structure &, const std::vector<LatticeSite> &) const;

  /// Returns a list of species associated with cluster space as chemical symbols.
  std::vector<std::vector<std::string>> getChemicalSymbols() const { return _chemicalSymbols; }

  /// Returns the cluster space size, i.e. the length of a cluster vector.
  size_t getClusterSpaceSize() { return _clusterVectorLength; }

  /// Returns the mapping between atomic numbers and the internal species enumeration scheme for each site.
  std::vector<std::unordered_map<int, int>> getSpeciesMaps() const { return _speciesMaps; }

  /// Returns the cluster product.
  double evaluateClusterProduct(const std::vector<int> &, const std::vector<int> &, const std::vector<int> &, const std::vector<int> &, const std::vector<int> &) const;

  /// Returns the default cluster function.
  double evaluateClusterFunction(const int, const int, const int) const;

  /// Computes permutations and multicomponent vectors of each orbit.
  void computeMultiComponentVectors();

  /// Removes orbits.
  void removeOrbits(std::vector<size_t> &);

  /// Merge orbits
  void mergeOrbits(const int index1, const int index2) { _primitiveOrbitList._orbits[index1] += _primitiveOrbitList._orbits[index2]; }

private:
  /// Primitive (prototype) structure.
  Structure _primitiveStructure;

  /// Primitive orbit list based on the structure and the global cutoffs
  OrbitList _primitiveOrbitList;

  /// Number of allowed components on each site of the primitive structure.
  std::vector<int> _numberOfAllowedSpeciesPerSite;

  /// Radial cutoffs by cluster order starting with pairs.
  std::vector<double> _clusterCutoffs;

  /// Species considered in this cluster space identified by atomic number.
  std::vector<int> _species;

  /// Map between atomic numbers and the internal species enumeration scheme for each site in the primitive structure.
  std::vector<std::unordered_map<int, int>> _speciesMaps;

  /// The allowed chemical symbols on each site in the primitive structure.
  std::vector<std::vector<std::string>> _chemicalSymbols;

  std::vector<std::vector<ClusterVectorElementInfo>> _clusterVectorElementInfoList;

  int _clusterVectorLength;
};
