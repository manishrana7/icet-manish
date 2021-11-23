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
    for (size_t i = 0; i < _orbitList.primitiveStructure().size(); i++)
    {
        Vector3d position_i = _orbitList.primitiveStructure().positions().row(i);
        LatticeSite lattice_site = _orbitList.primitiveStructure().findLatticeSiteByPosition(position_i, _fractionalPositionTolerance);
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
    _primitiveToSupercellMap.clear();

    std::set<Vector3d, Vector3dCompare> uniqueCellOffsets;

    // Map all sites
    for (size_t i = 0; i < _supercell->size(); i++)
    {
        Vector3d position_i = _supercell->positions().row(i);

        LatticeSite primitive_site = _orbitList.primitiveStructure().findLatticeSiteByPosition(position_i, _fractionalPositionTolerance);

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

    if (_uniquePrimcellOffsets.size() != _supercell->size() / _orbitList.primitiveStructure().size())
    {
        std::ostringstream msg;
        msg << "Wrong number of unitcell offsets found (LocalOrbitListGenerator::mapSitesAndFindCellOffsets)." << std::endl;
        msg << "Expected: " << _supercell->size() / _orbitList.primitiveStructure().size() << std::endl;
        msg << "Found:    " << _uniquePrimcellOffsets.size();
        throw std::runtime_error(msg.str());
    }
    std::sort(_uniquePrimcellOffsets.begin(), _uniquePrimcellOffsets.end(), Vector3dCompare());
}

/**
@details Generates and returns the local orbit list with the input index.
@param offset Cell offset in multiples of primitive lattice vectors
@param selfContained
    If this orbit list will be used on its own to calculate local cluster vectors or
    differences in cluster vector, this parameter needs to be true.
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
    std::shared_ptr<Structure> primitiveStructure = std::make_shared<Structure>(_orbitList.primitiveStructure());
    OrbitList localOrbitList = OrbitList(_orbitList.primitiveStructure());
    for (const Orbit &orbit : _orbitList.orbits())
    {
        // Copy the orbit.
        Orbit supercellOrbit = orbit;

        // If this orbit list will be used standalone for calculating
        // local cluster vectors or cluster vector differences,
        // we need to add clusters that include the present cell offset,
        // but would otherwise belong to the local orbit list of
        // another cell offset.
        if (selfContained)
        {
            // We will loop over the clusters in the orbit and add more
            // clusters inside the loop, so we first extract the original
            // clusters to avoid modifying the list we are looping over.
            std::vector<Cluster> clusters = supercellOrbit.clusters();
            for (auto cluster : clusters)
            {
                // Extract all versions of the clusters for which the
                // original cluster has been translated such that one
                // of the sites sits in the {0, 0, 0} cell offset.
                std::vector<std::vector<LatticeSite>> translatedSiteGroups = _orbitList.getSitesTranslatedToUnitcell(cluster.latticeSites(), false);
                for (auto translatedSites : translatedSiteGroups)
                {
                    // Only add clusters that are not duplicates of previus clusters.
                    // false or true here does not seem to matter.
                    if (!supercellOrbit.contains(translatedSites, true))
                    {
                        supercellOrbit.addCluster(Cluster(translatedSites, primitiveStructure));
                    }
                }
            }
        }

        // Translate all clusters of the new orbit.
        supercellOrbit.translate(offset);

        // Technically we should use the fractional position tolerance
        // corresponding to the cell metric of the supercell structure.
        // This is, however, not uniquely defined. Moreover, the difference
        // would only matter for very large supercells. We (@angqvist,
        // @erikfransson, @erhart) therefore decide to defer this issue
        // until someone encounters the problem in a practical situation.
        // In principle, one should not handle coordinates (floats) at this
        // level anymore. Rather one should transform any (supercell)
        // structure into an effective representation in terms of lattice
        // sites before any further operations.
        supercellOrbit.transformToSupercell(_supercell, _primitiveToSupercellMap, _fractionalPositionTolerance);
        localOrbitList.addOrbit(supercellOrbit);
    }
    return localOrbitList;
}

/// Generates the complete orbit list (the sum of all local orbit lists).
OrbitList LocalOrbitListGenerator::getFullOrbitList()
{
    OrbitList orbitList = OrbitList(*_supercell);
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
