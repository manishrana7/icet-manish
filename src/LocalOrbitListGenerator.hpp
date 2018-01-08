#pragma once
#include <vector>
#include <unordered_map>
#include "LatticeSite.hpp"
#include "Structure.hpp"
#include "OrbitList.hpp"
#include "Vector3dCompare.hpp"

/**

This is a small class that has a:

orbit list (from primitive structure)
supercell
list of unique primitive cell offsets that the supercell span
the primToSupercellMap


you can query this object with

///Give number of possible local orbit lists one can make from current supercell and the primitive orbit list
size_t number_of_possible_local_orbit_lists();

///Generate the orbit list from the primitive offset with count i
OrbitList getLocalOrbitList(int i);

std::vector<Vector3d> getUniqueOffsets() const;

std::vector<Vector3d> primToSupercellMap() const;

///clears primToSupercellMap and unique offsets
void reset();

etc...
*/

class LocalOrbitListGenerator
{
  public:
    LocalOrbitListGenerator(const OrbitList &, const Structure &);

    ///generate and returns the local orbit list with the input index
    OrbitList generateLocalOrbitList(const unsigned int ) ;

    ///generate and returns the local orbit list with the input offset (require that the offset is in uniquecell offset?)
    OrbitList generateLocalOrbitList(const Vector3d & ) ;

    /// Generate the full orbit list from this structure
    OrbitList generateFullOrbitList();

    //clears the unordered_map and the vector
    void clear();

    ///Returns the number of unique offsets
    size_t getUniqueOffsetsCount() const
    {
        return _uniquePrimcellOffsets.size();
    }

    ///Return the primitive lattice neighbor to supercell latticeneigbhor map
    std::unordered_map<LatticeSite, LatticeSite> getPrimToSupercellMap() const
    {
        return _primToSupercellMap;
    }

    ///Returns the unique primitive cells
    std::vector<Vector3d> getUniquePrimcellOffsets() const
    {
        return _uniquePrimcellOffsets;
    }



  private:

   /**
    Maps supercell positions to reference to the primitive cell and find unique primitive cell offsets
    Will loop through all sites in supercell and map them to the primitive structures cell
    and find the unique primitive cell offsets
    */
    void mapSitesAndFindCellOffsets();

    /// Find the sub permutation matrices that maps the basis atoms onto the supercell
    void findPermutationMatrices();


    ///Primitive orbit list
    OrbitList _orbit_list;

    ///supercell structure from which the local orbit list will be based upon
    Structure _supercell;

    ///this maps a latticeNeighbor from the primitive and get the equivalent in supercell
    std::unordered_map<LatticeSite, LatticeSite> _primToSupercellMap;

    ///The unique offsets of the primitive cell required to "cover" the supercell
    std::vector<Vector3d> _uniquePrimcellOffsets;

    /// The sub permutation matrices that will together map the basis atoms unto the supercell.
    std::vector<Matrix3i> _subPermutationMatrices;
};
