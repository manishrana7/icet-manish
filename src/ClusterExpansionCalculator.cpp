#include "ClusterExpansionCalculator.hpp"

ClusterExpansionCalculator::ClusterExpansionCalculator(const ClusterSpace &clusterSpace, const Structure &structure)
{
    _clusterSpace = clusterSpace;
    _clusterSpace._orbitList.sort();
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
            // for(const auto latticeSites : orbit.getPermutedEquivalentSites())
            {
                eqSiteIndex++;

                std::vector<LatticeSite> primitiveEquivalentSites;
                for (const auto site : latticeSites)
                {
                    Vector3d sitePosition = _superCell.getPosition(site);
                    auto primitiveSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(sitePosition);
                    primitiveEquivalentSites.push_back(primitiveSite);
                }
                //No translated sites (1)
                // std::vector<std::vector<LatticeSite>> latticeSitesTranslated;
                // latticeSitesTranslated.push_back(primitiveEquivalentSites);

                // //use translated sites (2)
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

    _fullPrimitiveOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
    int orbitIndex = -1;
    for (auto orbit : orbitVector)
    {
        orbitIndex++;
        // orbit.setEquivalentSitesPermutations(permutations[orbitIndex]);
        _fullPrimitiveOrbitList.addOrbit(orbit);
    }
    _fullPrimitiveOrbitList.addPermutationInformationToOrbits(_clusterSpace.getOrbitList().getCol1(),
                                                              _clusterSpace.getOrbitList().getPermutationMatrix());
    // _fullPrimitiveOrbitList.sort();
    // std::cout<<"Full list"<<std::endl;
    // _fullPrimitiveOrbitList.print();
    // validateBasisAtomOrbitLists();
    _primToSupercellMap.clear();
    for (int i = 0; i < structure.size(); i++)
    {
        Vector3d localPosition = structure.getPositions().row(i);
        LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition);
        Vector3d offsetVector = localSite.unitcellOffset();

        if (_localOrbitlists.find(offsetVector) == _localOrbitlists.end())
        {
            _localOrbitlists[offsetVector] = _fullPrimitiveOrbitList.getLocalOrbitList(structure, offsetVector, _primToSupercellMap);
            // for (auto &orbit : _localOrbitlists[offsetVector]._orbitList)
            // {
            //     auto permutedSites = orbit.getPermutedEquivalentSites();
            //     // orbit._equivalentSites = permutedSites;
            // }
        }
    }
    // testRemovingSites();
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
    bool orderIntact = true;  // dont sort the clusters
    bool permuteSites = true; // count the clusters in the order they lie in equivalent sites

    ClusterCounts clusterCounts = ClusterCounts();

    Vector3d localPosition = _superCell.getPositions().row(index);

    LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition);

    // Calc offset for this site

    Vector3d positionVector = _clusterSpace.getPrimitiveStructure().getPositions().row(0) - _clusterSpace.getPrimitiveStructure().getPositions().row(localSite.index());
    Vector3d localIndexZeroPos = localPosition + positionVector;
    auto indexZeroLatticeSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localIndexZeroPos);

    Vector3d offsetVector = indexZeroLatticeSite.unitcellOffset();
    Vector3d offsetVectorLocal = localSite.unitcellOffset();

    OrbitList translatedOrbitList = _localOrbitlists[offsetVector];
    if (offsetVector[0] != offsetVectorLocal[0] || offsetVector[1] != offsetVectorLocal[1] || offsetVector[2] != offsetVectorLocal[2])
    {
        translatedOrbitList += _localOrbitlists[offsetVectorLocal];
        std::cout << "adding extra OL" << std::endl;
    }
    for (const auto orbit : translatedOrbitList.getOrbitList())
    {
        for (const auto sites : orbit.getEquivalentSites())
        {
            if (std::none_of(sites.begin(), sites.end(), [](const LatticeSite &ls) { return ls.unitcellOffset().norm() < 1e-4; }))
            {
                // std::cout<<"No zeroth offsets sites in translated orbit"<<std::endl;
                // std::cout<<"Orbit order: "<<orbit.getRepresentativeCluster().order()<<std::endl;
            }
        }
    }
    translatedOrbitList.removeSitesNotContainingIndex(index);    
    // Purge the orbitlist of all sites containing the ignored indices
    for (auto ignoredIndex : ignoredIndices)
    {
        translatedOrbitList.removeSitesContainingIndex(ignoredIndex, -1);
    }

    int orbitIndex = -1;
    std::cout<<"calculating index "<<index<<std::endl;
    std::cout<<"ignored indices: ";
    for(const auto ignIndex : ignoredIndices){std::cout<< ignIndex<<" ";}
    std::cout<<std::endl;
    for(const auto orbit : _localOrbitlists[offsetVectorLocal]._orbitList)
    {
        orbitIndex++;
        for(const auto sites : orbit.getEquivalentSites())
        {
            if(!translatedOrbitList._orbitList[orbitIndex].contains(sites,true))
            {
                std::cout<<"Pruned these sites"<<std::endl;
                for(const auto site: sites)
                {
                    site.print();
                }
                std::cout<<std::endl;
            }
        }
    }
    translatedOrbitList.print();

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
            // std::cout << "Allowed occupations: ";
            // for (auto a : allowedOccupations)
            // {
            // std::cout << a << " ";
            // }
            // std::cout << std::endl;
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

        // auto mcVectors = _clusterSpace.getOrbitList()._orbitList[i].getMultiComponentVectors(allowedOccupations);
        auto mcVectors = _fullPrimitiveOrbitList._orbitList[i].getMultiComponentVectors(allowedOccupations);

        // auto allowedPermutationsSet = _clusterSpace.getOrbitList()._orbitList[i].getAllowedSitesPermutations();
        auto allowedPermutationsSet = _fullPrimitiveOrbitList._orbitList[i].getAllowedSitesPermutations();

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
                std::cout << "Didnt find cluster" << std::endl;
                continue;
            }
            for (const auto &elementsCountPair : clusterMap.at(repCluster))
            {
                for (const auto &perm : elementPermutations[currentMCVectorIndex])
                {
                    auto permutedMCVector = icet::getPermutedVector(mcVector, perm);
                    // std::cout << "mcVector: ";
                    // for (auto a : mcVector)
                    // {
                    //     std::cout << a << " ";
                    // }
                    // std::cout << std::endl;

                    // std::cout << "permutedMCVector: ";
                    // for (auto a : permutedMCVector)
                    // {
                    //     std::cout << a << " ";
                    // }
                    // std::cout << std::endl;
                    // std::cout << "elementsCountPair.first: ";
                    // for (auto a : elementsCountPair.first)
                    // {
                    //     std::cout << a << " ";
                    // }
                    // std::cout << std::endl;
                    auto permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);

                    // std::cout << "permutedAllowedOccupations: ";
                    // for (auto a : permutedAllowedOccupations)
                    // {
                    //     std::cout << a << " ";
                    // }
                    // std::cout << std::endl;
                    // std::cout << " cluster product " << _clusterSpace.evaluateClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) << std::endl;

                    clusterVectorElement += _clusterSpace.evaluateClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            //multiplicity != speciesPermutations[currentMultiComponentVectorIndex].size()*_orbitList.getOrbit(i).getEquivalentSites().size()*structure.size()/_primitiveStructure.size()
            double realMultiplicity = (double)elementPermutations[currentMCVectorIndex].size() * (double)_clusterSpace._orbitList._orbitList[i]._equivalentSites.size() / (double)_clusterSpace.getPrimitiveStructure().size();
            clusterVectorElement /= ((double)realMultiplicity * (double)_superCell.size());
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
            for (int j = i + 1; j < orbit.getEquivalentSites().size(); j++)
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

///Test that the removing of sites work as expected
void ClusterExpansionCalculator::testRemovingSites()
{
    // first get OL for the ith index
    for (int i = 0; i < _superCell.size(); i++)
    {
        // int i = 0;
        OrbitList orbitList_i = getLocalOrbitList(i);
        orbitList_i.removeSitesNotContainingIndex(i);

        // then get the Ol for the jth index
        for (int j = i + 1; j < _superCell.size(); j++)
        {
            // int j = 1;
            OrbitList orbitList_j = getLocalOrbitList(j);
            orbitList_j.removeSitesNotContainingIndex(j);
            orbitList_j.removeSitesContainingIndex(i, -1);

            // Now make sure there are no clusters in ol_i that contain site index i or j that exist in ol_j
            // Or the opposite

            for (int k = 0; k < orbitList_i.size(); k++)
            {
                for (const auto sites : orbitList_i.getOrbit(k).getEquivalentSites())
                {
                    // if (std::any_of(sites.begin(), sites.end(), [=](const LatticeSite &ls) { return ls.index() == i; }))
                    // {
                    if (orbitList_j.getOrbit(k).contains(sites, true))
                    {
                        std::cout << " double counting in i,j OL. i,j= " << i << ", " << j << std::endl;
                        for (const auto site : sites)
                        {
                            site.print();
                        }
                        std::cout << "<<end" << std::endl;
                    }
                    // }
                }
            }
        }
    }
}

OrbitList ClusterExpansionCalculator::getLocalOrbitList(int index)
{
    Vector3d localPosition = _superCell.getPositions().row(index);
    LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition);
    Vector3d positionVector = _clusterSpace.getPrimitiveStructure().getPositions().row(0) - _clusterSpace.getPrimitiveStructure().getPositions().row(localSite.index());
    Vector3d localIndexZeroPos = localPosition + positionVector;
    auto indexZeroLatticeSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localIndexZeroPos);
    Vector3d offsetVector = indexZeroLatticeSite.unitcellOffset();
    OrbitList orbitList_i = _localOrbitlists[offsetVector];
    return orbitList_i;
}