#include "LocalOrbitListGenerator.hpp"

LocalOrbitListGenerator::LocalOrbitListGenerator(const OrbitList &primitiveOrbitList, const Structure &superCell) : _orbit_list(primitiveOrbitList), _supercell(superCell)
{
    _positionClosestToOrigin = getClosestToOrigin();
    mapSitesAndFindCellOffsets();
    
    // generateSmartOffsets();
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
        Vector3d primitive_position = _orbit_list.getPrimitiveStructure().getPositions().row(primitive_site.index());
        // Basically only append offsets to indices that correspond to the atom in the origin
        if( (primitive_position - _positionClosestToOrigin).norm() < 1e-5)
        {
            uniqueCellOffsets.insert(primitive_site.unitcellOffset());
        }
    }

    _uniquePrimcellOffsets.clear();

    _uniquePrimcellOffsets.assign(uniqueCellOffsets.begin(), uniqueCellOffsets.end());

    std::sort(_uniquePrimcellOffsets.begin(), _uniquePrimcellOffsets.end(), Vector3dCompare());
}



void LocalOrbitListGenerator::generateSmartOffsets()
{
    int expected_number_of_mappings = _supercell.size() / _orbit_list.getPrimitiveStructure().size();
    // std::cout<<"supercell size "<<_supercell.size()<<std::endl;
    // std::cout<<"primitive size "<< _orbit_list.getPrimitiveStructure().size()<<std::endl;
    // std::cout<<"Expected number of mappings "<<expected_number_of_mappings<<std::endl;
    // Check that primitive size is a multiple of supercell size
    if (fabs(expected_number_of_mappings - _supercell.size() / _orbit_list.getPrimitiveStructure().size()) > 0.01)
    {
        std::string errorMessage = "Can not translage cell with ";
        errorMessage += std::to_string(_orbit_list.getPrimitiveStructure().size()) + " ";
        errorMessage += "to a super cell with ";
        errorMessage += std::to_string(_supercell.size()) + " atoms";
        throw std::runtime_error(errorMessage);
    }

    // Check if easy solution exists
    if (_uniquePrimcellOffsets.size() == expected_number_of_mappings)
    {
        /// Translate the normal way
    }

    // Get all the mappings that matches a supercell position by offsetting each primitive atom by the unique offsets
    std::vector<LocalOrbitListGenerator::Prim2SuperMap> allIndividualMappings = getIndividualMappings();

    if (allIndividualMappings.size() != _supercell.size())
    {
        throw std::runtime_error("The number individual mappings were not correct for this supercell");
    }

    // Create the smart offsets
    std::vector<Prim2SuperMappings> smartOffsets; 
    for (int j = 0; j < expected_number_of_mappings; j++)
    {
        // std::cout<<"j: "<<j<<" expected_number_of_mappings:"<<expected_number_of_mappings<<std::endl;
        Prim2SuperMappings currentMappings = Prim2SuperMappings();
        for (int i = allIndividualMappings.size()-1; i >=0; i--)
        {
            // std::cout<<"i: "<< i <<" allIndividualMappings.size():"<<allIndividualMappings.size()<<std::endl;
            if (isCompatibleNewMapping(currentMappings, allIndividualMappings[i]))
            {
                currentMappings.push_back(allIndividualMappings[i]);
                allIndividualMappings.erase(allIndividualMappings.begin()+i);
            }
            if(currentMappings.size() == expected_number_of_mappings)
            {
                break;
            }
        }
        if(currentMappings.size() != 0)
        {
            smartOffsets.push_back(currentMappings);
        }
    }
    
    std::cout<<" "<<std::endl;
    for(int i=0; i < smartOffsets.size(); i++)
    {
        std::cout<<smartOffsets[i].size()<<" ("<<_orbit_list.getPrimitiveStructure().size()<<")"<< " ("<<i<<" /"<<smartOffsets.size()<<")"<<std::endl;
    }
    for(const auto mapping : smartOffsets)
    {
        if( mapping.size() != _orbit_list.getPrimitiveStructure().size())
        {
            std::string errorMessage = "The different set of mappings were the wrong size ";
            errorMessage += std::to_string(_orbit_list.getPrimitiveStructure().size()) + " != ";
            errorMessage += std::to_string(mapping.size());
            throw std::runtime_error(errorMessage);
        }
    }
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
    for (int i = 0; i < getUniqueOffsetsCount(); i++)
    {
        orbitList += generateLocalOrbitList(i);
    }
    return orbitList;
}

/// Get the vector of all the individual primtive  to super via offset
std::vector<LocalOrbitListGenerator::Prim2SuperMap> LocalOrbitListGenerator::getIndividualMappings() const
{
    std::vector<LocalOrbitListGenerator::Prim2SuperMap> individualMappings;
    for (const auto &offset : _uniquePrimcellOffsets)
    {
        std::vector<Vector3d> positions = icet::getOffsetPositions(_orbit_list.getPrimitiveStructure(), offset);
        for (int i = 0; i < positions.size(); i++)
        {
            std::vector<int> matchedIndices = findMatchingSupercellPositions(positions[i]);

            for (const auto superIndex : matchedIndices)
            {
                Prim2SuperMap p2sm = Prim2SuperMap(i, superIndex, offset);
                individualMappings.push_back(p2sm);
            }
        }
    }
    return individualMappings;
}

/// Find the indices of the supercell that are within 1e-3 of the position argument
std::vector<int> LocalOrbitListGenerator::findMatchingSupercellPositions(const Vector3d &position) const
{
    std::vector<int> matchedIndices;
    for (int i = 0; i < _supercell.getPositions().rows(); i++)
    {
        if ((position.transpose() - _supercell.getPositions().row(i)).norm() < 1e-3)
        {
            matchedIndices.push_back(i);
        }
    }
    return matchedIndices;
}

/// Test if next mapping is compatible to be added to curentMappings
bool LocalOrbitListGenerator::isCompatibleNewMapping(LocalOrbitListGenerator::Prim2SuperMappings currentMappings, LocalOrbitListGenerator::Prim2SuperMap mapping) const
{
    if (currentMappings.primitiveIndices.size() == 0)
    {
        return true;
    }
    if (currentMappings.isPrimitiveIndexInside(mapping.primitiveIndex))
    {
        return false;
    }
    if (currentMappings.isSuperIndexInside(mapping.superIndex))
    {
        return false;
    }


    /**
     *  Check that the distance between each primitve index to the candiate
     *  primitive index is the same as the
     * 
     * 
     */
    for(int i = 0; i < currentMappings.primitiveIndices.size(); i++)
    {
        int primitiveIndex = currentMappings.primitiveIndices[i];
        Vector3d offset = currentMappings.offsets[i];

        int superIndex = currentMappings.superIndices[i];
        Vector3d zeroOffset = {0., 0., 0.};

        double distNoOffset = _orbit_list.getPrimitiveStructure().getDistance(primitiveIndex,
                                                                           mapping.primitiveIndex,
                                                                           zeroOffset, zeroOffset);
        // auto superLatticeSite = _supercell.findLatticeSiteByPosition(pos1);
        // auto superLatticeSiteTrial = _supercell.findLatticeSiteByPosition(pos2);
        double distOffset = _orbit_list.getPrimitiveStructure().getDistance(primitiveIndex,
                                                                           mapping.primitiveIndex,
                                                                           offset, mapping.offset);
        double distPrim = _orbit_list.getPrimitiveStructure().getDistance(primitiveIndex,
                                                                           mapping.primitiveIndex,
                                                                           zeroOffset, zeroOffset);
        double distSuper = _supercell.getDistance(superIndex,mapping.superIndex, zeroOffset, zeroOffset);
        if (fabs(distNoOffset - distSuper) > 1e-2)
        {
            return false;
        }
    }
    return true;
}

//clears the unordered_map and the vector
void LocalOrbitListGenerator::clear()
{
    _primToSupercellMap.clear();
    _uniquePrimcellOffsets.clear();
}
