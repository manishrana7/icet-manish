#include "LocalOrbitListGenerator.hpp"

/**
@param orbitList orbit list for the underlying primitive cell
@param supercell supercell structure for which to set up the local orbit list generation
@param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates
*/
LocalOrbitListGenerator::LocalOrbitListGenerator(const OrbitList &orbitList,
                                                 std::shared_ptr<Structure> supercell,
                                                 const double fractionalPositionTolerance)
    : _orbitList(orbitList),
      _supercell(supercell),
      _fractionalPositionTolerance(fractionalPositionTolerance)
{
    _indexToClosestAtom = getIndexOfAtomClosestToOrigin();
    mapSitesAndFindCellOffsets();
}

/// @details This position is used for extracting unit cell offsets later on.
int LocalOrbitListGenerator::getIndexOfAtomClosestToOrigin()
{
    double distanceToOrigin = 1e6;
    int indexToClosestAtom;
    for (size_t i = 0; i < _orbitList.getPrimitiveStructure().size(); i++)
    {
        Vector3d position_i = _orbitList.getPrimitiveStructure().getPositions().row(i);
        LatticeSite lattice_site = _orbitList.getPrimitiveStructure().findLatticeSiteByPosition(position_i, _fractionalPositionTolerance);
        // @todo Can this be removed?
        if (lattice_site.unitcellOffset().norm() > FLOATTYPE_EPSILON)
        {
            continue;
        }
        if (position_i.norm() < distanceToOrigin)
        {
            distanceToOrigin = position_i.norm();
            indexToClosestAtom = i;
        }
    }
    return indexToClosestAtom;
}

/**
@details Maps supercell positions to reference to the primitive cell and find
unique primitive cell offsets. Loops through all sites in supercell and
map them to the primitive structures cell and find the unique primitive cell
offsets.
*/
void LocalOrbitListGenerator::mapSitesAndFindCellOffsets()
{
    _primToSupercellMap.clear();

    std::set<Vector3d, Vector3dCompare> uniqueCellOffsets;

    // Map all sites
    for (size_t i = 0; i < _supercell->size(); i++)
    {
        Vector3d position_i = _supercell->getPositions().row(i);

        LatticeSite primitive_site = _orbitList.getPrimitiveStructure().findLatticeSiteByPosition(position_i, _fractionalPositionTolerance);

        // @todo Can we just use zero and remove
        // the getIndexOfAtomClosestToOrigin function?
        if (primitive_site.index() == _indexToClosestAtom)
        {
            uniqueCellOffsets.insert(primitive_site.unitcellOffset());
        }
    }

    // If empty: add zero offset
    if (uniqueCellOffsets.size() == 0)
    {
        Vector3d zeroVector = {0.0, 0.0, 0.0};
        uniqueCellOffsets.insert(zeroVector);
    }

    _uniquePrimcellOffsets.clear();

    _uniquePrimcellOffsets.assign(uniqueCellOffsets.begin(), uniqueCellOffsets.end());

    if (_uniquePrimcellOffsets.size() != _supercell->size() / _orbitList.getPrimitiveStructure().size())
    {
        std::ostringstream msg;
        msg << "Wrong number of unitcell offsets found (LocalOrbitListGenerator::mapSitesAndFindCellOffsets)." << std::endl;
        msg << "Expected: " << _supercell->size() / _orbitList.getPrimitiveStructure().size() << std::endl;
        msg << "Found:    " << _uniquePrimcellOffsets.size();
        throw std::runtime_error(msg.str());
    }
    std::sort(_uniquePrimcellOffsets.begin(), _uniquePrimcellOffsets.end(), Vector3dCompare());
}

/**
@details Generates and returns the local orbit list with the input index.
@param CXXXXXXX
*/
OrbitList LocalOrbitListGenerator::getLocalOrbitList(const Vector3d &offset, bool selfContained = false)
{
    if (std::find(_uniquePrimcellOffsets.begin(), _uniquePrimcellOffsets.end(),
                  offset) == _uniquePrimcellOffsets.end())
    {
        std::ostringstream msg;
        msg << "The offset " << offset << "was not found in _uniquePrimcellOffsets(LocalOrbitListGenerator::getLocalOrbitList)" << std::endl;
        throw std::out_of_range(msg.str());
    }
    return _orbitList.getLocalOrbitList(_supercell, offset, _primToSupercellMap, _fractionalPositionTolerance, selfContained);
}

/// Generates the complete orbit list (the sum of all local orbit lists).
OrbitList LocalOrbitListGenerator::getFullOrbitList()
{
    OrbitList orbitList = OrbitList();
    for (auto offset : _uniquePrimcellOffsets)
    {
        orbitList += getLocalOrbitList(offset);
    }

    if (orbitList.size() != _orbitList.size())
    {
        std::ostringstream msg;
        msg << "Full orbitlist size is not the same as local orbitlist size (LocalOrbitListGenerator::getFullOrbitList)" << std::endl;
        msg << " full orbitlist size: " << orbitList.size() << std::endl;
        msg << " local orbitlist size: " << _orbitList.size() << std::endl;
        throw std::runtime_error(msg.str());
    }
    return orbitList;
}
