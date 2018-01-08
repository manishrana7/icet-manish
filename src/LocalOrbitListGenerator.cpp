#include "LocalOrbitListGenerator.hpp"

LocalOrbitListGenerator::LocalOrbitListGenerator(const OrbitList &primitiveOrbitList, const Structure &superCell) : _orbit_list(primitiveOrbitList), _supercell(superCell)
{
    mapSitesAndFindCellOffsets();
}




/**
    Maps supercell positions to reference to the primitive cell and find unique primitive cell offsets
    Will loop through all sites in supercell and map them to the primitive structures cell
    and find the unique primitive cell offsets
    */
void LocalOrbitListGenerator::mapSitesAndFindCellOffsets()
{
    _primToSupercellMap.clear();

    std::set<Vector3d, Vector3dCompare> uniqueCellOffsets;

    //map all sites
    for (size_t i = 0; i < _supercell.size(); i++)
    {
        Vector3d position_i = _supercell.getPositions().row(i);

        LatticeSite primitive_site = _orbit_list.getPrimitiveStructure().findLatticeSiteByPosition(position_i);
        // LatticeSite super_site = _supercell.findLatticeSiteByPosition(position_i);

        // _primToSupercellMap[primitive_site] = super_site;
        uniqueCellOffsets.insert(primitive_site.unitcellOffset());
    }

    _uniquePrimcellOffsets.clear();

    _uniquePrimcellOffsets.assign(uniqueCellOffsets.begin(), uniqueCellOffsets.end());

    std::sort(_uniquePrimcellOffsets.begin(), _uniquePrimcellOffsets.end(), Vector3dCompare());
}

///generate and returns the local orbit list with the input index
OrbitList LocalOrbitListGenerator::generateLocalOrbitList(const unsigned int index)
{
    if (index >= _uniquePrimcellOffsets.size())
    {
        std::string errMsg = "Error: attempting to generateLocalOrbitList with index " + std::to_string(index) + " when size of unique offsets are: " + std::to_string(_uniquePrimcellOffsets.size()) + ".";
        throw std::out_of_range(errMsg);
    }

    return _orbit_list.getLocalOrbitList(_supercell, _uniquePrimcellOffsets[index], _primToSupercellMap);
}

///generate and returns the local orbit list with the input offset (require that the offset is in uniquecell offset?)
OrbitList LocalOrbitListGenerator::generateLocalOrbitList(const Vector3d &primOffset)
{
    auto find = std::find(_uniquePrimcellOffsets.begin(), _uniquePrimcellOffsets.end(), primOffset);
    if (find == _uniquePrimcellOffsets.end())
    {
        std::cout << "Warning: generating local orbit list with offset not found in _uniquePrimcellOffsets" << std::endl;
    }

    return _orbit_list.getLocalOrbitList(_supercell, primOffset, _primToSupercellMap);
}


/// Generate the complete orbit list (the sum of all local orbit lists)
OrbitList LocalOrbitListGenerator::generateFullOrbitList()
{
    OrbitList orbitList = OrbitList();
    for(int i = 0; i < getUniqueOffsetsCount(); i++)
    {
        orbitList += generateLocalOrbitList(i);
    }
    return orbitList;


}


//clears the unordered_map and the vector
void LocalOrbitListGenerator::clear()
{
    _primToSupercellMap.clear();
    _uniquePrimcellOffsets.clear();
}
