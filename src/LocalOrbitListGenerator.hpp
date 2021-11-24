#pragma once

#include <iomanip>
#include <unordered_map>
#include <vector>

#include "LatticeSite.hpp"
#include "OrbitList.hpp"
#include "Structure.hpp"
#include "VectorOperations.hpp"

/**

This is a small class that has a:

orbit list (from primitive structure)
supercell
list of unique primitive cell offsets that the supercell span
the primitiveToSupercellMap


you can query this object with

///Generate the orbit list from the primitive offset with count i
OrbitList getLocalOrbitList(int i);

std::vector<Vector3d> getUniqueOffsets() const;

std::vector<Vector3d> primitiveToSupercellMap() const;

etc...
*/

class LocalOrbitListGenerator
{
public:
    /// Constructor.
    LocalOrbitListGenerator(const OrbitList &, std::shared_ptr<Structure>, const double);

    /// Generates and returns the local orbit list with the input index.
    OrbitList getLocalOrbitList(const Vector3d &, bool);

    /// Generates the full orbit list from this structure.
    OrbitList getFullOrbitList();

    /// Returns the number of unique offsets.
    size_t getNumberOfUniqueOffsets() const { return _uniquePrimcellOffsets.size(); }

    /// Returns the primitive lattice neighbor to supercell lattice neigbhor map.
    std::unordered_map<LatticeSite, LatticeSite> getMapFromPrimitiveToSupercell() const { return _primitiveToSupercellMap; }

    /// Returns the unique primitive cells
    std::vector<Vector3d> getUniquePrimitiveCellOffsets() const { return _uniquePrimcellOffsets; }

private:
    /// Maps supercell positions to reference.
    void mapSitesAndFindCellOffsets();

    /// Primitive orbit list.
    OrbitList _primitiveOrbitList;

    /// Primitive structure on which the supercell is based.
    std::shared_ptr<Structure> _primitiveStructure;

    /// Supercell structure on which the local orbit list will be based.
    std::shared_ptr<Structure> _supercell;

    /// Maps a lattice site from the primitive cell to an equivalent lattice site in the supercell.
    std::unordered_map<LatticeSite, LatticeSite> _primitiveToSupercellMap;

    /// Finds the index of the atom that is closest to the origin.
    int getIndexOfAtomClosestToOrigin();

    /// @todo Add description
    size_t _indexToClosestAtom;

    /// The unique offsets of the primitive cell required to "cover" the supercell.
    std::vector<Vector3d> _uniquePrimcellOffsets;

    /// Tolerance applied when comparing positions in Cartesian coordinates.
    double _fractionalPositionTolerance;
};
