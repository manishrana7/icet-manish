#include "ClusterExpansionCalculator.hpp"

ClusterExpansionCalculator::ClusterExpansionCalculator(const ClusterSpace &clusterSpace, const Structure &structure)
{
    _clusterSpace = clusterSpace;    
    _superCell = structure;
    
    LocalOrbitListGenerator _theLog = LocalOrbitListGenerator(clusterSpace.getOrbitList(), _superCell);
    size_t uniqueOffsets = _theLog.getNumberOfUniqueOffsets();
    int numberOfOrbits = _clusterSpace._orbitList.size();
    std::vector<Orbit> orbitVector;
    for (const auto orbit : clusterSpace._orbitList._orbitList)
    {
        orbitVector.push_back(Orbit(orbit.getRepresentativeCluster()));
    }

    // Permutations for the sites in the orbits
    std::vector<std::vector<std::vector<int>>> permutations(numberOfOrbits);

    /* Strategy to construct the "full" primitive orbitlists
    
    We first fill up a std::vector<Orbit> orbitVector
    where vector<orbit> is essentially an orbitlist.
    
    The existing methods to construct the full orbitlist is to loop over all the
    unique offsets.

    We loop over each local orbitlist (by looping over offsetIndex)
    The local orbitlist is retrieved here:
        `_theLog.getLocalOrbitList(offsetIndex).getOrbitList()`

    Then for each group of latticesites in orbit.equivalentSites() we add them
    to orbitVector[i][orbitIndex] if the latticesites has a site with offset [0, 0, 0].

    */

    for (int offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
    {
        int orbitIndex = -1;
        // This orbit is a local orbit related to the supercell
        for (const auto orbit : _theLog.getLocalOrbitList(offsetIndex).getOrbitList())
        {
            orbitIndex++;

            auto orbitPermutations = orbit.getEquivalentSitesPermutations();

            int eqSiteIndex = -1;

            for (const auto latticeSites : orbit.getEquivalentSites())
            {
                eqSiteIndex++;

                std::vector<LatticeSite> primitiveEquivalentSites;
                for (const auto site : latticeSites)
                {
                    Vector3d sitePosition = _superCell.getPosition(site);
                    auto primitiveSite = _clusterSpace._primitiveStructure.findLatticeSiteByPosition(sitePosition);
                    primitiveEquivalentSites.push_back(primitiveSite);
                }
                std::vector<std::vector<LatticeSite>> latticeSitesTranslated = _clusterSpace._orbitList.getSitesTranslatedToUnitcell(primitiveEquivalentSites, false);

                for (auto latticesitesPrimTrans : latticeSitesTranslated)
                {
                    if (std::any_of(latticesitesPrimTrans.begin(), latticesitesPrimTrans.end(), [=](LatticeSite ls) { return (ls.unitcellOffset()).norm() < 1e-4; }))
                    {
                        // false or true here seems to not matter
                        if (!orbitVector[orbitIndex].contains(latticesitesPrimTrans, true))
                        {
                            orbitVector[orbitIndex].addEquivalentSites(latticesitesPrimTrans);
                            permutations[orbitIndex].push_back(orbitPermutations[eqSiteIndex]);
                        }
                    }
                }
            }
        }
    }
    /// Now create the full primitive orbitlist using the vector<orbit>
    _fullPrimitiveOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
    int orbitIndex = -1;
    for (auto orbit : orbitVector)
    {
        orbitIndex++;
        _fullPrimitiveOrbitList.addOrbit(orbit);
    }
    _fullPrimitiveOrbitList.addPermutationInformationToOrbits(_clusterSpace.getOrbitList().getCol1(),
                                                              _clusterSpace.getOrbitList().getPermutationMatrix());

    _primToSupercellMap.clear();
    _indexToOffset.clear();
    for (int i = 0; i < structure.size(); i++)
    {
        Vector3d localPosition = structure.getPositions().row(i);
        LatticeSite localSite = _clusterSpace._primitiveStructure.findLatticeSiteByPosition(localPosition);
        Vector3d offsetVector = localSite.unitcellOffset();
        _indexToOffset[i] = offsetVector;

        if (_localOrbitlists.find(offsetVector) == _localOrbitlists.end())
        {
            _localOrbitlists[offsetVector] = _fullPrimitiveOrbitList.getLocalOrbitList(structure, offsetVector, _primToSupercellMap);
            for (auto &orbit : _localOrbitlists[offsetVector]._orbitList)
            {
                auto permutedSites = orbit.getPermutedEquivalentSites();
                orbit._equivalentSites = permutedSites;
            }
        }
    }
}

/**
@details This constructs a cluster vector that only considers clusters that contain the input index.
@param occupations the occupation vector for the supercell
@param index the local index of the supercell 
@param ignoredIndices a vector of indices which have already had their local energy calculated. This is required to input so that no double counting occurs.
*/

std::vector<double> ClusterExpansionCalculator::getLocalClusterVector(const std::vector<int> &occupations, int index, std::vector<int> ignoredIndices)
{
    _superCell.setAtomicNumbers(occupations);

    if (occupations.size() != _superCell.size())
    {
        throw std::runtime_error("Input occupations and internal supercell structure mismatch in size");
    }

    for (auto ignoreIndex : ignoredIndices)
    {
        if (ignoreIndex >= _superCell.size())
        {
            throw std::runtime_error("index larger than Input structure size in method ClusterExpansionCalculator::getLocalClusterVector");
        }
    }
    
    // dont sort the clusters
    bool orderIntact = true;   

    // count the clusters in the order they lie in equivalent sites
    // since these sites are already in the permuted order
    bool permuteSites = false; 

    // Remove all that doesn't contain index regardless of offset?
    bool removeGhostIndexNotContain = true;

    // Remove all ignored indices regardless of offset?
    bool removeGhostIndexContain = false;

    ClusterCounts clusterCounts = ClusterCounts();

    // Get one of the translated orbitlists
    OrbitList translatedOrbitList = _localOrbitlists[_indexToOffset[index]];

    // Remove sites not containing the local index
    if(_clusterSpace._primitiveStructure.size()>1)
    {
        translatedOrbitList.removeSitesNotContainingIndex(index, removeGhostIndexNotContain);
    }

    // Purge the orbitlist of all sites containing the ignored indices
    for (auto ignoredIndex : ignoredIndices)
    {
        translatedOrbitList.removeSitesContainingIndex(ignoredIndex, removeGhostIndexContain);
    }

    // Count clusters and get cluster count map
    clusterCounts.countOrbitList(_superCell, translatedOrbitList, orderIntact, permuteSites);
    const auto clusterMap = clusterCounts.getClusterCounts();

    // Finally begin occupying the cluster vector
    int orbitIndex = -1;
    std::vector<double> clusterVector;
    clusterVector.push_back(1.0 / _superCell.size());
    for (size_t i = 0; i < _fullPrimitiveOrbitList.size(); i++)
    {
        Cluster repCluster = _fullPrimitiveOrbitList._orbitList[i]._representativeCluster;
        std::vector<int> allowedOccupations;

        if (i >= _clusterSpace._orbitList.size())
        {
            std::cout << _fullPrimitiveOrbitList.size() << " >= " << _clusterSpace._orbitList.size() << std::endl;
            throw std::runtime_error("index i larger than cs.orbit_list.size() in ClusterExpansionCalculator::getLocalClusterVector");
        }
        try
        {
            allowedOccupations = _clusterSpace.getNumberOfAllowedSpeciesBySite(_clusterSpace._primitiveStructure, _clusterSpace._orbitList._orbitList[i].getRepresentativeSites());
        }
        catch (const std::exception &e)
        {
            std::cout << e.what() << std::endl;
            throw std::runtime_error("Failed getting allowed occupations in genereteClusterVector");
        }

        // Skip rest if any sites aren't active sites (i.e. allowed occupation < 2)
        if (std::any_of(allowedOccupations.begin(), allowedOccupations.end(), [](int allowedOccupation) { return allowedOccupation < 2; }))
        {
            continue;
        }

        const auto &mcVectors = _clusterSpace._multiComponentVectors[i];
        const auto &elementPermutations = _clusterSpace._sitePermutations[i];
        repCluster.setTag(i);
        for (int currentMCVectorIndex = 0; currentMCVectorIndex < _clusterSpace._multiComponentVectors[i].size(); currentMCVectorIndex++)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            auto clusterFind = clusterMap.find(repCluster);
            if (clusterFind == clusterMap.end())
            {
                clusterVector.push_back(0);
                continue;
            }
            for (const auto &elementsCountPair : clusterMap.at(repCluster))
            {
                for (const auto &perm : _clusterSpace._sitePermutations[i][currentMCVectorIndex])
                {
                    const auto &permutedMCVector = icet::getPermutedVector(mcVectors[currentMCVectorIndex], perm);
                    const auto &permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);

                    clusterVectorElement += _clusterSpace.evaluateClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }            
            double realMultiplicity = (double)_clusterSpace._sitePermutations[i][currentMCVectorIndex].size() * (double)_clusterSpace._orbitList._orbitList[i]._equivalentSites.size() / (double)_clusterSpace._primitiveStructure.size();
            clusterVectorElement /= ((double)realMultiplicity * (double)_superCell.size());
            clusterVector.push_back(clusterVectorElement);
        }
    }
    return clusterVector;
}
