#pragma once

#include <vector>

#include "Cluster.hpp"
#include "LatticeSite.hpp"
#include "ManyBodyNeighborList.hpp"
#include "Orbit.hpp"
#include "Structure.hpp"
#include "Symmetry.hpp"
#include "VectorOperations.hpp"
#include "VectorOperations.hpp"

/**
This class serves as a container for a sorted list of orbits and provides associated functionality.
*/

class OrbitList
{
public:
    /// Empty constructor.
    OrbitList(){};

    /// Constructs orbit list from a set of neighbor lists, a matrix of equivalent sites, and a structure.
    OrbitList(const Structure &,
              const std::vector<std::vector<LatticeSite>> &,
              const std::vector<std::vector<std::vector<LatticeSite>>> &,
              const double);

    /// Sort orbit list.
    void sort(const double);

    /// Adds an orbit.
    void addOrbit(const Orbit &orbit);

    /// Adding-to-existing (+=) operator.
    OrbitList &operator+=(const OrbitList &);

    /// Returns the orbit of the given index.
    const Orbit &getOrbit(unsigned int) const;

    /// Returns the local orbit list for a site.
    OrbitList getLocalOrbitList(std::shared_ptr<Structure>,
                                const Vector3d &,
                                std::unordered_map<LatticeSite, LatticeSite> &,
                                const double) const;

    // @todo Add description.
    void addColumnsFromMatrixOfEquivalentSites(std::vector<std::vector<std::vector<LatticeSite>>> &,
                                               std::unordered_set<std::vector<int>, VectorHash> &,
                                               const std::vector<int> &,
                                               bool) const;

    /// Clears the orbit list.
    void clear() { _orbits.clear(); }

    /// Returns the number of orbits.
    size_t size() const { return _orbits.size(); }

    // Returns the first column of the matrix of equivalent sites.
    std::vector<LatticeSite> getReferenceLatticeSites(bool = true) const;

    // Returns rows of the matrix of equivalent sites that match the lattice sites.
    std::vector<int> getIndicesOfEquivalentLatticeSites(const std::vector<LatticeSite> &,
                                                        bool = true) const;

    // Returns true if the cluster includes at least on site from the unit cell at the origin.
    bool validCluster(const std::vector<LatticeSite> &) const;

    // @todo Add description.
    std::vector<LatticeSite> translateSites(const std::vector<LatticeSite> &,
                                            const unsigned int) const;

    /// @todo Add description.
    std::vector<std::vector<LatticeSite>> getSitesTranslatedToUnitcell(const std::vector<LatticeSite> &,
                                                                       bool sortit = true) const;

    /// @todo Add description.
    std::vector<std::pair<std::vector<LatticeSite>, std::vector<int>>> getMatchesInMatrixOfEquivalenSites(const std::vector<std::vector<LatticeSite>> &) const;

    /// @todo Add description.
    void transformSiteToSupercell(LatticeSite &,
                                  const Structure &,
                                  std::unordered_map<LatticeSite, LatticeSite> &,
                                  const double) const;

    /// Returns all columns from the given rows in matrix of equivalent sites
    std::vector<std::vector<LatticeSite>> getAllColumnsFromRow(const std::vector<int> &,
                                                               bool,
                                                               bool) const;

    /// @todo Clarify description.
    /// First construct rows_sort = sorted(rows)  then returns true/false if rows_sort exists in taken_rows
    bool isRowsTaken(const std::unordered_set<std::vector<int>,
                                              VectorHash> &,
                     std::vector<int>) const;

    /// Finds and returns sites in first column of matrix of equivalent sites along with their unit cell translated indistinguishable sites.
    std::vector<std::vector<LatticeSite>> getAllColumnsFromCluster(const std::vector<LatticeSite> &) const;

    /// Returns the first column of the matrix of equivalent sites used to construct the orbit list.
    std::vector<LatticeSite> getFirstColumnOfMatrixOfEquivalentSites() const { return _referenceLatticeSites; }

    /// Returns the matrix of equivalent sites used to construct the orbit list.
    std::vector<std::vector<LatticeSite>> getMatrixOfEquivalentSites() const { return _matrixOfEquivalentSites; }

    /// Removes an orbit identified by its index.
    void removeOrbit(const size_t);

    /// Removes all orbits that have inactive sites.
    void removeInactiveOrbits(const Structure &);

    /// Returns the orbits in this orbit list.
    const std::vector<Orbit> getOrbits() const { return _orbits; }

    /// Returns the primitive structure.
    Structure getPrimitiveStructure() const { return _primitiveStructure; }

    /// Sets primitive structure.
    void setPrimitiveStructure(const Structure &primitive) { _primitiveStructure = primitive; }

    /// Merge two orbits.
    void mergeOrbits(int index1, int index2);

private:
    /// Contains all the orbits in the orbit list.
    std::vector<Orbit> _orbits;

    /// @todo Add description.
    std::vector<LatticeSite> _referenceLatticeSites;

    /// Matrix of equivalent sites.
    std::vector<std::vector<LatticeSite>> _matrixOfEquivalentSites;

    /// Primitive structure for which orbit list was constructed.
    Structure _primitiveStructure;

    /// Create an orbit by deducing the proper permutations
    Orbit createOrbit(const std::vector<std::vector<LatticeSite>> &);

    std::shared_ptr<Structure> _primitiveStructurePtr;
};
