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
    /// This collects the information that primtiveIndex is mapped to superIndex via the offset
    struct Prim2SuperMap
    {
        Prim2SuperMap(int index, int superIndex, Vector3d offset) : primitiveIndex(index), superIndex(superIndex), offset(offset)
        {
        }
        int primitiveIndex;
        int superIndex;
        Vector3d offset;
    };

    /// This contains one set of mappings that take the entire primitive to one part of the supercell
    struct Prim2SuperMappings
    {
        std::vector<int> primitiveIndices;
        std::vector<int> superIndices;
        std::vector<Vector3d> offsets;

        Prim2SuperMap getMap(int index)
        {
            if (index >= primitiveIndices.size())
            {
                throw std::out_of_range("Index out of range in Prim2SuperMap getMap(int index)");
            }
            Prim2SuperMap p2sm = Prim2SuperMap(primitiveIndices[index], superIndices[index], offsets[index]);
            return p2sm;
        }
        size_t size() const
        {
            return primitiveIndices.size();
        }
        void push_back(const Prim2SuperMap p2sm)
        {
            primitiveIndices.push_back(p2sm.primitiveIndex);
            superIndices.push_back(p2sm.superIndex);
            offsets.push_back(p2sm.offset);
        }
        bool isPrimitiveIndexInside(int index)
        {
            auto find = std::find(primitiveIndices.begin(), primitiveIndices.end(), index);
            return find != primitiveIndices.end();
        }
        bool isSuperIndexInside(int index)
        {
            auto find = std::find(superIndices.begin(), superIndices.end(), index);
            return find != superIndices.end();
        }

    };

  
    LocalOrbitListGenerator(const OrbitList &, const Structure &);

    ///generate and returns the local orbit list with the input index
    OrbitList generateLocalOrbitList(const unsigned int);

    ///generate and returns the local orbit list with the input offset (require that the offset is in uniquecell offset?)
    OrbitList generateLocalOrbitList(const Vector3d &);

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

    /// Generate smart offsets @ TODO change name
    void generateSmartOffsets();

    ///Primitive orbit list
    OrbitList _orbit_list;

    ///supercell structure from which the local orbit list will be based upon
    Structure _supercell;

    ///this maps a latticeNeighbor from the primitive and get the equivalent in supercell
    std::unordered_map<LatticeSite, LatticeSite> _primToSupercellMap;

    /// Find the position of the atom that is closest to the origin.
    /// This position is used to extracting unit cell offsets later on
    Vector3d getClosestToOrigin() const
    {
       Vector3d closestToOrigin;
       double distanceToOrigin = 1e6;
       for(int i=0; i < _orbit_list.getPrimitiveStructure().size(); i++)
       {
           Vector3d position_i = _orbit_list.getPrimitiveStructure().getPositions().row(i);
        if( (position_i.norm()) < distanceToOrigin)
        {
            distanceToOrigin = position_i.norm();
            closestToOrigin = position_i;
        }
       }
       return closestToOrigin;
    }

    ///
    Vector3d _positionClosestToOrigin;
    ///The unique offsets of the primitive cell required to "cover" the supercell
    std::vector<Vector3d> _uniquePrimcellOffsets;

    /// The sub permutation matrices that will together map the basis atoms unto the supercell.
    std::vector<Matrix3i> _subPermutationMatrices;

    /// Get the vector of all the individual primtive  to super via offset
    std::vector<Prim2SuperMap> getIndividualMappings() const;

    /// Find the indices of the supercell that are within 1e-3 of the position argument
    std::vector<int> findMatchingSupercellPositions(const Vector3d &position) const;

    /// Test if next mapping is compatible to be added to curentMappings
    bool isCompatibleNewMapping(LocalOrbitListGenerator::Prim2SuperMappings currentMappings, LocalOrbitListGenerator::Prim2SuperMap mapping) const;


};
