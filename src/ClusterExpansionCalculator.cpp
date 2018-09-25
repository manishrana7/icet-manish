#include "ClusterExpansionCalculator.hpp"

ClusterExpansionCalculator::ClusterExpansionCalculator(const ClusterSpace &clusterSpace, const Structure &structure)
{
    _clusterSpace = clusterSpace;
    _superCell = structure;
    // Set up stuff.
    _theLog = LocalOrbitListGenerator(clusterSpace.getOrbitList(), _superCell);
    size_t uniqueOffsets = _theLog.getNumberOfUniqueOffsets();
    int numberOfOrbits = _clusterSpace._orbitList.size();
    Vector3d zeroOffset = {0.0, 0.0, 0.0};
    int primitiveSize = _clusterSpace.getPrimitiveStructure().size();
    std::vector<Orbit> orbitVector;
    // std::cout<<"Primitve"<<std::endl;
    // clusterSpace.getOrbitList().print();
    for (const auto orbit : clusterSpace._orbitList._orbitList)
    {
        orbitVector.push_back(Orbit(orbit.getRepresentativeCluster()));
    }

    // Permutations for the sites in the orbits
    std::vector<std::vector<std::vector<int>>> permutations(numberOfOrbits);

    /* Strategy to construct the "full" primitive orbitlists
    
    We first fill up a std::vector<std::vector<Orbit>> orbitVector
    where the inner vector<orbit> is essentially an orbitlist.
    the outer vector is over basis atoms. When the entire vector<vector<orbit>> 
    is constructed we create a vector<orbitlists>

    The existing methods to construct the full orbitlist is to loop over all the
    unique offsets.

    We loop over each local orbitlist (by looping over offsetIndex)
    The local orbitlist is retrieved here:
        `_theLog.getLocalOrbitList(offsetIndex).getOrbitList()`

    Then for each group of latticesites in orbit.equivalentSites() we add them
    to orbitVector[i][orbitIndex] if the latticesites has the basis atom `i` in them.

    */

   // Map lattice indices to other lattice indices that are equivalent
   std::map<int, int> singletMap;


    for (int offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
    {
        int orbitIndex = -1;
        // This orbit is a local orbit related to the supercell
        for (const auto orbit : _theLog.getLocalOrbitList(offsetIndex).getOrbitList())
        {
            orbitIndex++;

            auto orbitPermutations = orbit.getEquivalentSitesPermutations();

            int eqSiteIndex = -1;

            for (const auto latticeSites : orbit._equivalentSites)
            {
                eqSiteIndex++;

                for (int i = 0; i < primitiveSize; i++)
                {
                    Vector3d primPos = _clusterSpace.getPrimitiveStructure().getPositions().row(i);
                    LatticeSite primitiveSite_i = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(primPos);
                    LatticeSite superCellEquivalent = _superCell.findLatticeSiteByPosition(primPos);

                    std::vector<LatticeSite> primitiveEquivalentSites;
                    for (const auto site : latticeSites)
                    {
                        Vector3d sitePosition = _superCell.getPosition(site);
                        auto primitiveSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(sitePosition);
                        primitiveEquivalentSites.push_back(primitiveSite);
                    }

                    std::vector<std::vector<LatticeSite>> latticeSitesTranslated = _clusterSpace._orbitList.getSitesTranslatedToUnitcell(primitiveEquivalentSites);                                        
                    for (auto latticesitesPrimTrans : latticeSitesTranslated)
                    {
                        if (std::any_of(latticesitesPrimTrans.begin(), latticesitesPrimTrans.end(), [=](LatticeSite ls) { return (ls.unitcellOffset()).norm() < 1e-4; }))
                        {
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
    }

    _fullPrimitiveOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
    int orbitIndex = -1;
    for (auto &orbit : orbitVector)
    {
        orbitIndex++;
        orbit.setEquivalentSitesPermutations(permutations[orbitIndex]);
        _fullPrimitiveOrbitList.addOrbit(orbit);
    }
    // std::cout<<"Full list"<<std::endl;
    // _fullPrimitiveOrbitList.print();
    validateBasisAtomOrbitLists();
    for (int i = 0; i < structure.size(); i++)
    {
        Vector3d localPosition = structure.getPositions().row(i);
        LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition);
        Vector3d offsetVector = localSite.unitcellOffset();

        if (_localOrbitlists.find(offsetVector) == _localOrbitlists.end())
        {
            _localOrbitlists[offsetVector] = _fullPrimitiveOrbitList.getLocalOrbitList(structure, offsetVector, _primToSupercellMap);
            for (auto &orbit : _localOrbitlists[offsetVector]._orbitList)
            {
                auto permutedSites = orbit.getPermutedEquivalentSites();
                // orbit._equivalentSites = permutedSites;
            }
        }
    }
}

/**
@details This constructs a cluster vector that considers only clusters with local to the input index.
@param numberOfAllowedSpecies number of allowed components for each site of the primitive structure
@param chemicalSymbols chemical symbol for each site
@param orbitList list of orbits for the primitive structure
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

    int dprint = 0;
    bool orderIntact = true;   // dont sort the clusters
    bool permuteSites = true; // count the clusters in the order they lie in equivalent sites

    ClusterCounts clusterCounts = ClusterCounts();

    Vector3d localPosition = _superCell.getPositions().row(index);

    LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition);

    // Calc offset for this site

    Vector3d positionVector = _clusterSpace.getPrimitiveStructure().getPositions().row(0)-_clusterSpace.getPrimitiveStructure().getPositions().row(localSite.index());
    Vector3d localIndexZeroPos = localPosition + positionVector;
    auto indexZeroLatticeSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localIndexZeroPos);

    Vector3d offsetVector = indexZeroLatticeSite.unitcellOffset();
    // Vector3d offsetVector = localSite.unitcellOffset();

    OrbitList translatedOrbitList = _localOrbitlists[offsetVector];
    translatedOrbitList.removeSitesNotContainingIndex(index);
    // Purge the orbitlist of all sites containing the ignored indices
    for (auto ignoredIndex : ignoredIndices)
    {
        translatedOrbitList.removeSitesContainingIndex(ignoredIndex);
    }

    clusterCounts.countOrbitList(_superCell, translatedOrbitList, orderIntact, permuteSites);
    const auto clusterMap = clusterCounts.getClusterCounts();

    std::vector<double> clusterVector;
    clusterVector.push_back(1.0 / _superCell.size());
    // Finally begin occupying the cluster vector
    for (size_t i = 0; i < _fullPrimitiveOrbitList.size(); i++)
    {
        Cluster repCluster = _fullPrimitiveOrbitList._orbitList[i].getRepresentativeCluster();
        std::vector<int> allowedOccupations;

        if (i >= _clusterSpace._orbitList.size())
        {
            std::cout << _fullPrimitiveOrbitList.size() << " >= " << _clusterSpace._orbitList.size() << std::endl;
            throw std::runtime_error("index i larger than cs.orbit_list.size() in ClusterExpansionCalculator::getLocalClusterVector");
        }
        try
        {
            allowedOccupations = _clusterSpace.getNumberOfAllowedSpeciesBySite(_clusterSpace.getPrimitiveStructure(), _clusterSpace._orbitList._orbitList[i].getRepresentativeSites());
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

        auto mcVectors = _clusterSpace._orbitList._orbitList[i].getMultiComponentVectors(allowedOccupations);
        auto allowedPermutationsSet = _clusterSpace._orbitList._orbitList[i].getAllowedSitesPermutations();
        auto elementPermutations = _clusterSpace.getMultiComponentVectorPermutations(mcVectors, i);
        repCluster.setTag(i);
        int currentMCVectorIndex = 0;
        for (const auto &mcVector : mcVectors)
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
                for (const auto &perm : elementPermutations[currentMCVectorIndex])
                {
                    auto permutedMCVector = icet::getPermutedVector(mcVector, perm);
                    auto permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);
                    clusterVectorElement += _clusterSpace.evaluateClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            int realMultiplicity = elementPermutations[currentMCVectorIndex].size()*_clusterSpace._orbitList._orbitList[i]._equivalentSites.size()/_clusterSpace.getPrimitiveStructure().size();
            clusterVectorElement /= ((double)realMultiplicity*(double)_superCell.size());
            clusterVector.push_back(clusterVectorElement);
            currentMCVectorIndex++;
        }
    }
    return clusterVector;
}

// std::vector<double> getLocalClusterVector(const Structure &, const std::vector<int>)

void ClusterExpansionCalculator::validateBasisAtomOrbitLists()
{

    for (const auto orbit : _fullPrimitiveOrbitList._orbitList)
    {
        for (int i = 0; i < orbit._equivalentSites.size(); i++)
        {
            for (int j = i + 1; j < orbit._equivalentSites.size(); j++)
            {
                auto sites_i = orbit.getEquivalentSites()[i];
                std::sort(sites_i.begin(), sites_i.end());
                auto sites_j = orbit.getEquivalentSites()[j];
                std::sort(sites_j.begin(), sites_j.end());
                if (std::equal(sites_i.begin(), sites_i.end(), sites_j.begin()))
                {
                    throw std::runtime_error("Two eq. sites in an orbit were equal.");
                }
            }
        }
    }
}