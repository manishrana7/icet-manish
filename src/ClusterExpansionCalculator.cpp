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
    std::vector<std::vector<Orbit>> orbitVector(uniqueOffsets, std::vector<Orbit>(numberOfOrbits));
    OrbitList basisOrbitList = OrbitList();

    for (int offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
    {
        int orbitIndex = 0;
        for (const auto orbit : _theLog.generateLocalOrbitList(offsetIndex).getOrbitList())
        {
            auto orbitPermutations = orbit.getEquivalentSitesPermutations();

            int eqSiteIndex = 0;
            for (const auto latticeSites : orbit.getEquivalentSites())
            {
                for (int i = 0; i < primitiveSize; i++)
                {
                    Vector3d primPos = _clusterSpace.getPrimitiveStructure().getPositions().row(i);
                    LatticeSite primitiveSite_i = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(primPos);
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
                    }

                    orbitIndex++;
                } //end orbit loop
            }
        }
    }

            // old code below

            for (int i = 0; i < _clusterSpace.getPrimitiveStructure().size(); i++)
            {
                OrbitList basisOrbitList = OrbitList();
                basisOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
                Vector3d primPos = _clusterSpace.getPrimitiveStructure().getPositions().row(i);
                LatticeSite primitiveSite_i = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(primPos);

                LatticeSite superCellEquivalent = _theLog.getPrimToSupercellMap()[primitiveSite_i];

                int orbitIndex = 0;
                std::vector<Orbit> orbitVector();
                for (int offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
                {

                    // for (const auto orbit : fullOrbitList.getOrbitList())
                    for (const auto orbit : _theLog.generateLocalOrbitList(offsetIndex).getOrbitList())
                    {
                        Orbit basisOrbit = Orbit(orbit.getRepresentativeCluster());
                        auto orbitPermutations = orbit.getEquivalentSitesPermutations();
                        std::vector<std::vector<int>> permutations;
                        int eqSiteIndex = 0;

                        for (const auto latticeSites : orbit.getEquivalentSites())
                        {

                            if (std::find(latticeSites.begin(), latticeSites.end(), superCellEquivalent) != latticeSites.end())
                            {
                                std::vector<LatticeSite> primitiveEquivalentSites;
                                for (const auto site : latticeSites)
                                {
                                    Vector3d sitePosition = _superCell.getPosition(site);
                                    auto primitiveSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(sitePosition);
                                    primitiveEquivalentSites.push_back(primitiveSite);
                                }
                                // basisOrbit.addEquivalentSites(latticeSites);
                                basisOrbit.addEquivalentSites(primitiveEquivalentSites);
                                if (orbit.getClusterSize() > 1)
                                {
                                    if (eqSiteIndex >= orbitPermutations.size())
                                    {
                                        throw std::runtime_error("eqSiteIndex and orbitPermutations are not the same size");
                                    }
                                    permutations.push_back(orbitPermutations[eqSiteIndex]);
                                }
                            }
                            eqSiteIndex++;
                        }
                        if (permutations.size() != basisOrbit.size() and orbit.getClusterSize() > 1)
                        {
                            throw std::runtime_error("Permurations and basisorbit are not the same size");
                        }

                        basisOrbit.setEquivalentSitesPermutations(permutations);
                        basisOrbitList.addOrbit(basisOrbit);
                    }
                }
                _basisAtomOrbitList.push_back(basisOrbitList);
            }
            // std::cout<<"Done in init"<<std::endl;
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
                auto repCluster = _basisAtomOrbitList[basisIndex].getOrbit(i).getRepresentativeCluster();
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

                        // throw std::runtime_error("represtentative cluster not found in clusterMap");
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
                    clusterVectorElement /= ((double)multiplicity * structure.size());
                    clusterVector.push_back(clusterVectorElement * mcVector.size());

                    currentMCVectorIndex++;
                }
            }
            return clusterVector;
        }
