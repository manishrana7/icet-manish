#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "LocalOrbitListGenerator.hpp"
#include "OrbitList.hpp"
#include "PeriodicTable.hpp"
#include "Structure.hpp"
#include "VectorOperations.hpp"

/**
@brief This struct keeps track of information pertaining to a specific element
       in the cluster vector.
*/
struct ClusterVectorElementInfo
{
  /// A multi-component vector contains the indices of the point functions
  /// (only non-trivial if the number of components are more than two)
  std::vector<int> multiComponentVector;

  /// Site permutations describe how the sites in the cluster can be re-ordered
  std::vector<std::vector<int>> sitePermutations;

  /// The index in the cluster vector for this element
  int clusterVectorIndex;

  /// Multiplicity for this cluster vector element
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
  ClusterSpace(std::vector<std::vector<std::string>> &, std::shared_ptr<OrbitList>, const double, const double);

  /// Returns the cluster vector corresponding to the input structure.
  std::vector<double> getClusterVector(const Structure &, const double) const;

  /// Returns information concerning the association between orbits and multi-component vectors.
  std::pair<int, std::vector<int>> getMultiComponentVectorsByOrbit(const unsigned int);

  /// Returns the entire primitive orbit list.
  const OrbitList &getPrimitiveOrbitList() const { return *_primitiveOrbitList; }

  /// Returns the multi-component (MC) vector permutations for each MC vector in the set of input vectors.
  std::vector<std::vector<std::vector<int>>> getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &, const int) const;

  /// Returns the cluster vector given the orbit list and a structure.
  const std::vector<double> getClusterVectorFromOrbitList(const OrbitList &, const std::shared_ptr<Structure>, const double firstElement = 1.0, const int flipIndex = -1, const int newOccupation = -1) const;

  /// Returns the cutoff for each order.
  std::vector<double> getCutoffs() const { return _clusterCutoffs; }

  /// Returns the primitive structure.
  const Structure &primitiveStructure() const { return _primitiveStructure; }

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

  /// Computes permutations and multi-component vectors of each orbit.
  void computeMultiComponentVectors();

  /// Removes orbits.
  void removeOrbits(std::vector<size_t> &);

  /// Merge orbits.
  void mergeOrbits(const int index1, const int index2) { _primitiveOrbitList->mergeOrbits(index1, index2); }

private:
  /// Primitive (prototype) structure.
  Structure _primitiveStructure;

  /// Primitive orbit list based on the structure and the cutoffs.
  std::shared_ptr<OrbitList> _primitiveOrbitList;

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

  /// Information about each cluster vector element (multiplicity, multi-component vectors etc.).
  std::vector<std::vector<ClusterVectorElementInfo>> _clusterVectorElementInfoList;

  /// Length of the cluster vector.
  int _clusterVectorLength;
};
