#include "ClusterExpansionCalculator.hpp"

ClusterExpansionCalculator::ClusterExpansionCalculator(const ClusterSpace &clusterSpace, const Structure &structure)
{
    _clusterSpace = clusterSpace;
    _superCell = structure;
    _basisAtomOrbitList.clear();
    // Set up stuff.
    _theLog = LocalOrbitListGenerator(clusterSpace.getOrbitList(), _superCell);
    size_t uniqueOffsets = _theLog.getUniqueOffsetsCount();
    int numberOfOrbits = _clusterSpace.getOrbitList().size();
    Vector3d zeroOffset = {0.0, 0.0, 0.0};
    int primitiveSize = _clusterSpace.getPrimitiveStructure().size();
    /// New code here
    std::vector<std::vector<Orbit>> orbitVector(primitiveSize);
    for (int i = 0; i < primitiveSize; i++)
    {
        for (const auto orbit : clusterSpace.getOrbitList().getOrbitList())
        {
            orbitVector[i].push_back(Orbit(orbit.getRepresentativeCluster()));
        }
    }

    std::cout<<"Primitive size: "<<primitiveSize<<std::endl;
    // Permutations for the sites in the orbits
    // vector<vector<int> > v(10, vector<int>(10,1));
    std::vector<std::vector<std::vector<std::vector<int>>>> permutations(primitiveSize, std::vector<std::vector<std::vector<int>>>(numberOfOrbits));
    OrbitList basisOrbitList = OrbitList();
    /* Strategy to construct the "full" primitive orbitlists
    
    We first full up a std::vector<std::vector<Orbit>> orbitVector
    where the inner vector<orbit> is essentially an orbitlist.
    the outer vector is over basis atoms. When the entire vector<vector<orbit>> 
    is constructed we create a vector<orbitlists> 

    The existing methods to construct the full orbitlist is to loop over all the
    unique offsets.

    We loop over each local orbitlist (by looping over offsetIndex)
    The local orbitlist is retrieved here:
        `_theLog.generateLocalOrbitList(offsetIndex).getOrbitList()`

    Then for each group of latticesites in orbit.equivalentSites() we add them
    to orbitVector[i][orbitIndex] if the latticesites has the basis atom `i` in them.

    */
    for (int offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
    {
        int orbitIndex = -1;
        for (const auto orbit : _theLog.generateLocalOrbitList(offsetIndex).getOrbitList())
        {
            orbitIndex++;

            auto orbitPermutations = orbit.getEquivalentSitesPermutations();

            int eqSiteIndex = -1;
            for (const auto latticeSites : orbit.getEquivalentSites())
            {
                eqSiteIndex++;
                for (int i = 0; i < primitiveSize; i++)
                {

                    Vector3d primPos = _clusterSpace.getPrimitiveStructure().getPositions().row(i);
                    LatticeSite primitiveSite_i = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(primPos);
                    LatticeSite superCellEquivalent = _theLog.getPrimToSupercellMap()[primitiveSite_i];

                    if (std::find(latticeSites.begin(), latticeSites.end(), superCellEquivalent) != latticeSites.end())
                    {
                        std::vector<LatticeSite> primitiveEquivalentSites;
                        for (const auto site : latticeSites)
                        {
                            Vector3d sitePosition = _superCell.getPosition(site);
                            auto primitiveSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(sitePosition);
                            primitiveEquivalentSites.push_back(primitiveSite);
                        }
                        orbitVector[i][orbitIndex].addEquivalentSites(primitiveEquivalentSites);
                        permutations[i][orbitIndex].push_back(orbitPermutations[eqSiteIndex]);
                    }
                }
            }
        }
    }

    for (int i = 0; i < primitiveSize; i++)
    {
        OrbitList orbitList = OrbitList();
        orbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
        int orbitIndex = -1;
        for (auto &orbit : orbitVector[i])
        {
            orbitIndex++;
            orbit.setEquivalentSitesPermutations(permutations[i][orbitIndex]);
            orbitList.addOrbit(orbit);
        }
        _basisAtomOrbitList.push_back(orbitList);
    }

    // // old code below

    // for (int i = 0; i < _clusterSpace.getPrimitiveStructure().size(); i++)
    // {
    //     OrbitList basisOrbitList = OrbitList();
    //     basisOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
    //     Vector3d primPos = _clusterSpace.getPrimitiveStructure().getPositions().row(i);
    //     LatticeSite primitiveSite_i = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(primPos);

    //     LatticeSite superCellEquivalent = _theLog.getPrimToSupercellMap()[primitiveSite_i];

    //     int orbitIndex = 0;
    //     std::vector<Orbit> orbitVector();
    //     for (int offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
    //     {

    //         // for (const auto orbit : fullOrbitList.getOrbitList())
    //         for (const auto orbit : _theLog.generateLocalOrbitList(offsetIndex).getOrbitList())
    //         {
    //             Orbit basisOrbit = Orbit(orbit.getRepresentativeCluster());
    //             auto orbitPermutations = orbit.getEquivalentSitesPermutations();
    //             std::vector<std::vector<int>> permutations;
    //             int eqSiteIndex = 0;

    //             for (const auto latticeSites : orbit.getEquivalentSites())
    //             {

    //                 if (std::find(latticeSites.begin(), latticeSites.end(), superCellEquivalent) != latticeSites.end())
    //                 {
    //                     std::vector<LatticeSite> primitiveEquivalentSites;
    //                     for (const auto site : latticeSites)
    //                     {
    //                         Vector3d sitePosition = _superCell.getPosition(site);
    //                         auto primitiveSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(sitePosition);
    //                         primitiveEquivalentSites.push_back(primitiveSite);
    //                     }
    //                     // basisOrbit.addEquivalentSites(latticeSites);
    //                     basisOrbit.addEquivalentSites(primitiveEquivalentSites);
    //                     if (orbit.getClusterSize() > 1)
    //                     {
    //                         if (eqSiteIndex >= orbitPermutations.size())
    //                         {
    //                             throw std::runtime_error("eqSiteIndex and orbitPermutations are not the same size");
    //                         }
    //                         permutations.push_back(orbitPermutations[eqSiteIndex]);
    //                     }
    //                 }
    //                 eqSiteIndex++;
    //             }
    //             if (permutations.size() != basisOrbit.size() and orbit.getClusterSize() > 1)
    //             {
    //                 throw std::runtime_error("Permurations and basisorbit are not the same size");
    //             }

    //             basisOrbit.setEquivalentSitesPermutations(permutations);
    //             basisOrbitList.addOrbit(basisOrbit);
    //         }
    //     }
    //     _basisAtomOrbitList.push_back(basisOrbitList);
    // }
    validateBasisAtomOrbitLists();
    std::cout << "Done in init" << std::endl;
}

std::vector<double> ClusterExpansionCalculator::getLocalClusterVector(const Structure &structure, const int index)
{

    if (structure.size() != _superCell.size())
    {
        throw std::runtime_error("Input structure and internal supercell structure mismatch in size");
    }

    if (index >= structure.size())
    {
        throw std::runtime_error("index larger than Input structure size in method ClusterExpansionCalculator::getLocalClusterVector");
    }

    int dprint = 0;
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster

    ClusterCounts clusterCounts = ClusterCounts();

    Vector3d localPosition = structure.getPositions().row(index);

    LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition);

    int basisIndex = localSite.index();
    if (basisIndex >= _basisAtomOrbitList.size())
    {
        throw std::runtime_error("basisIndex and _basisAtomOrbitList are not the same size");
    }
    // auto primToSupercellMap = _theLog.getPrimToSupercellMap();

    for (const auto orbit : _basisAtomOrbitList[basisIndex].getOrbitList())
    {
        for (const auto eqSites : orbit.getEquivalentSites())
        {
            for (const auto site : eqSites)
            {
                if (site.index() > _basisAtomOrbitList[basisIndex].getPrimitiveStructure().size())
                {
                    throw std::runtime_error("lattice site index in orbit is larger than prim.size()");
                }
                if (site.index() > _clusterSpace.getPrimitiveStructure().size())
                {
                    throw std::runtime_error("lattice site index in orbit is larger than prim.size()");
                }
            }
        }
    }

    auto translatedOrbitList = _basisAtomOrbitList[basisIndex].getLocalOrbitList(structure, localSite.unitcellOffset(), _primToSupercellMap);
    clusterCounts.countOrbitList(structure, translatedOrbitList, orderIntact);

    // for( const auto orbit : _basisAtomOrbitList[basisIndex].getOrbitList() )
    // {
    //     std::cout<<"orbit order "<< orbit.getClusterSize()<< " orbit size "<<orbit.size()<<std::endl;
    // }

    const auto clusterMap = clusterCounts.getClusterCounts();
    std::vector<double> clusterVector;
    clusterVector.push_back(1.0 / structure.size());
    // Finally begin occupying the cluster vector
    for (size_t i = 0; i < _basisAtomOrbitList[basisIndex].size(); i++)
    {
        // std::cout<<"i="<<i<<std::endl;
        Cluster repCluster = _basisAtomOrbitList[basisIndex].getOrbit(i).getRepresentativeCluster();
        std::vector<int> allowedOccupations;

        if (i >= _clusterSpace.getOrbitList().size())
        {
            std::cout << _basisAtomOrbitList[basisIndex].size() << " >= " << _clusterSpace.getOrbitList().size() << std::endl;
            throw std::runtime_error("index i larger than cs.orbit_list.size() in ClusterExpansionCalculator::getLocalClusterVector");
        }
        try
        {
            allowedOccupations = _clusterSpace.getAllowedOccupations(_clusterSpace.getPrimitiveStructure(), _clusterSpace.getOrbitList().getOrbit(i).getRepresentativeSites());
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

        auto mcVectors = _clusterSpace.getOrbitList().getOrbit(i).getMCVectors(allowedOccupations);
        auto allowedPermutationsSet = _clusterSpace.getOrbitList().getOrbit(i).getAllowedSitesPermutations();
        auto elementPermutations = _clusterSpace.getMCVectorPermutations(mcVectors, i);
        repCluster.setClusterTag(i);
        int currentMCVectorIndex = 0;
        for (const auto &mcVector : mcVectors)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            // std::cout<<"do clusterMap.at(repCluster) "<<std::endl;
            // repCluster.print();

            auto clusterFind = clusterMap.find(repCluster);

            if (clusterFind == clusterMap.end())
            {
                // std::cout<<"repcluster: ";
                // repCluster.print();
                // std::cout<<"cluster tag: "<<repCluster.tag()<<std::endl;
                // std::cout<<"clusterMap size "<< clusterMap.size()<<std::endl;
                // std::cout<<"orbitlist size "<< _basisAtomOrbitList[basisIndex].size()<<std::endl;

                // std::cout<<" cluster map cluster key tags: "<<std::endl;
                // for(const auto mapPair : clusterMap)
                // {
                //     std::cout<<mapPair.first.tag()<< " ";
                // }
                // std::cout<<std::endl;

                throw std::runtime_error("represtentative cluster not found in clusterMap");
                clusterVector.push_back(0);
                continue;
            }
            for (const auto &elementsCountPair : clusterMap.at(repCluster))
            {
                // std::cout<<"done clusterMap.at(repCluster) "<<std::endl;

                // TODO check if allowedOccupations should be permuted as well.
                for (const auto &perm : elementPermutations[currentMCVectorIndex])
                {
                    auto permutedMCVector = icet::getPermutedVector(mcVector, perm);
                    auto permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);
                    clusterVectorElement += _clusterSpace.getClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            int realMultiplicity = _clusterSpace.getOrbitList().getOrbitList()[i].getEquivalentSites().size();
            clusterVectorElement /= ((double)structure.size() * realMultiplicity  );
            clusterVector.push_back(clusterVectorElement * repCluster.order());

            currentMCVectorIndex++;
        }
    }
    return clusterVector;
}

void ClusterExpansionCalculator::validateBasisAtomOrbitLists()
{

    for (const auto orbitList : _basisAtomOrbitList)
    {
        for (const auto orbit : orbitList.getOrbitList())
        {
            for (int i = 0; i < orbit.getEquivalentSites().size(); i++)
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
}