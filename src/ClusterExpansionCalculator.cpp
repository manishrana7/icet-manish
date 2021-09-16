#include "ClusterExpansionCalculator.hpp"

ClusterExpansionCalculator::ClusterExpansionCalculator(const ClusterSpace &clusterSpace,
                                                       const Structure &structure,
                                                       const double fractionalPositionTolerance)
{
    _clusterSpace = clusterSpace;
    _supercell = structure;
    LocalOrbitListGenerator LOLG = LocalOrbitListGenerator(clusterSpace.getOrbitList(), _supercell, fractionalPositionTolerance);
    size_t uniqueOffsets = LOLG.getNumberOfUniqueOffsets();
    int numberOfOrbits = _clusterSpace._orbitList.size();
    std::vector<Orbit> orbitVector;
    _fullOrbitList = LOLG.getFullOrbitList();
    for (const auto orbit : clusterSpace._orbitList._orbits)
    {
        orbitVector.push_back(Orbit(orbit.getRepresentativeCluster()));
    }

    // Permutations for the clusters in the orbits
    std::vector<std::vector<std::vector<int>>> permutations(numberOfOrbits);

    /* Strategy for constructing the "full" primitive orbit lists.

    First we fill up a `std::vector<Orbit> orbitVector`,
    where `vector<orbit>` is essentially an orbit list.

    The existing method for constructing the _full_ orbit list proceeds
    by looping over all local orbit lists with `LocalOrbitListGenerator` and
    adding the sites to the local orbit list.

    Now we do something similar by looping over each local orbit list
    (by looping over `offsetIndex`).
    The local orbitlist is retrieved here:
        `LOLG.getLocalOrbitList(offsetIndex).getOrbits()`

    Then for each orbit `orbitIndex` in `LOLG.getLocalOrbitList(offsetIndex).getOrbits()`
    each group of lattice sites in `orbit.equivalentSites()` is added to
    `orbitVector[orbitIndex]` if the lattice sites have a site with offset `[0, 0, 0]`.

    When the full primitive orbit list is used to create a local orbit list for
    site `index` in the supercell it should thus contain all lattice sites that
    contain `index`.
    */

    for (size_t offsetIndex = 0; offsetIndex < uniqueOffsets; offsetIndex++)
    {
        int orbitIndex = -1;
        // This orbit is a local orbit related to the supercell
        for (const auto orbit : LOLG.getLocalOrbitList(offsetIndex).getOrbits())
        {
            orbitIndex++;

            auto orbitPermutations = orbit.getPermutationsOfEquivalentClusters();

            int eqSiteIndex = -1;

            for (const auto cluster : orbit.getEquivalentClusters())
            {
                eqSiteIndex++;

                std::vector<LatticeSite> primitiveEquivalentSites;
                for (const auto site : cluster)
                {
                    Vector3d sitePosition = _supercell.getPosition(site);
                    auto primitiveSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(sitePosition, fractionalPositionTolerance);
                    primitiveEquivalentSites.push_back(primitiveSite);
                }
                std::vector<std::vector<LatticeSite>> translatedClusters = _clusterSpace._orbitList.getSitesTranslatedToUnitcell(primitiveEquivalentSites, false);

                for (auto translatedCluster : translatedClusters)
                {
                    if (std::any_of(translatedCluster.begin(), translatedCluster.end(), [=](LatticeSite ls)
                                    { return (ls.unitcellOffset()).norm() < fractionalPositionTolerance; }))
                    {
                        // false or true here does not seem to matter
                        if (!orbitVector[orbitIndex].contains(translatedCluster, true))
                        {
                            orbitVector[orbitIndex].addEquivalentCluster(translatedCluster);
                            permutations[orbitIndex].push_back(orbitPermutations[eqSiteIndex]);
                        }
                    }
                }
            }
        }
    }

    // Now create the full primitive orbit list using the vector<orbit>
    _fullPrimitiveOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());
    int orbitIndex = -1;
    for (auto orbit : orbitVector)
    {
        orbitIndex++;
        orbit.setPermutationsOfEquivalentClusters(permutations[orbitIndex]);
        _fullPrimitiveOrbitList.addOrbit(orbit);
    }

    // Calculate the permutation for each orbit in this orbit list.
    // This is normally done in the constructor but since we made one manually
    // we have to do this ourself.
    // _fullPrimitiveOrbitList.addPermutationInformationToOrbits(_clusterSpace.getOrbitList().getFirstColumnOfMatrixOfEquivalentSites(),
    //                                                           _clusterSpace.getOrbitList().getMatrixOfEquivalentSites());

    _primToSupercellMap.clear();
    _indexToOffset.clear();

    // Precompute all possible local orbitlists for this supercell and map it to the offset
    for (size_t i = 0; i < structure.size(); i++)
    {
        Vector3d localPosition = structure.getPositions().row(i);
        LatticeSite localSite = _clusterSpace.getPrimitiveStructure().findLatticeSiteByPosition(localPosition, fractionalPositionTolerance);
        Vector3d offsetVector = localSite.unitcellOffset();
        _indexToOffset[i] = offsetVector;

        if (_localOrbitlists.find(offsetVector) == _localOrbitlists.end())
        {
            _localOrbitlists[offsetVector] = _fullPrimitiveOrbitList.getLocalOrbitList(structure, offsetVector, _primToSupercellMap, fractionalPositionTolerance);

            // Set equivalent cluster equal to the permuted clusters so no permutation is required in the orbit list counting.
            for (auto &orbit : _localOrbitlists[offsetVector]._orbits)
            {
                auto permutedClusters = orbit.getPermutedEquivalentClusters();
                orbit._equivalentClusters = permutedClusters;
            }
        }
    }
}

/**
@details Occupy cluster vector based on information currently in _clusterCounts
@param First element of the cluster vector
*/
std::vector<double> ClusterExpansionCalculator::_occupyClusterVector(const double firstElement)
{
    std::vector<double> clusterVector;
    clusterVector.push_back(firstElement);
    for (size_t i = 0; i < _fullPrimitiveOrbitList.size(); i++)
    {
        Cluster representativeCluster = _fullPrimitiveOrbitList._orbits[i]._representativeCluster;
        auto representativeSites = _clusterSpace._orbitList._orbits[i].getSitesOfRepresentativeCluster();
        std::vector<int> allowedOccupations;

        if (i >= _clusterSpace._orbitList.size())
        {
            std::cout << _fullPrimitiveOrbitList.size() << " >= " << _clusterSpace._orbitList.size() << std::endl;
            throw std::runtime_error("Index i larger than cs.orbit_list.size() (ClusterExpansionCalculator::getLocalClusterVector)");
        }
        try
        {
            allowedOccupations = _clusterSpace.getNumberOfAllowedSpeciesBySite(_clusterSpace._primitiveStructure, representativeSites);
        }
        catch (const std::exception &e)
        {
            std::cout << e.what() << std::endl;
            throw std::runtime_error("Failed getting allowed occupations (ClusterExpansionCalculator::getLocalClusterVector)");
        }

        // Skip the rest if any of the sites are inactive (i.e. allowed occupation < 2)
        if (std::any_of(allowedOccupations.begin(), allowedOccupations.end(), [](int allowedOccupation)
                        { return allowedOccupation < 2; }))
        {
            continue;
        }

        std::vector<int> indicesOfRepresentativeSites;
        for (const auto site : representativeSites)
        {
            indicesOfRepresentativeSites.push_back(site.index());
        }

        const auto &mcVectors = _clusterSpace._multiComponentVectors[i];
        representativeCluster.setTag(i);

        /// Loop over all multi component vectors for this orbit
        for (size_t currentMCVectorIndex = 0; currentMCVectorIndex < _clusterSpace._multiComponentVectors[i].size(); currentMCVectorIndex++)
        {
            double clusterVectorElement = 0;

            auto clusterFind = _clusterCounts._clusterCounts.find(representativeCluster);

            /// Push back zero if nothing was counted for this orbit
            if (clusterFind == _clusterCounts._clusterCounts.end())
            {
                clusterVector.push_back(0);
                continue;
            }

            std::vector<int> permutedMCVector;
            std::vector<int> permutedAllowedOccupations;
            std::vector<int> permutedRepresentativeIndices;
            /// Loop over all the counts for this orbit
            for (const auto &elementsCountPair : _clusterCounts._clusterCounts.at(representativeCluster))
            {
                /// Loop over all equivalent permutations for this orbit and mc vector
                for (const auto &perm : _clusterSpace._sitePermutations[i][currentMCVectorIndex])
                {
                    /// Permute the mc vector and the allowed occupations
                    permutedMCVector = icet::getPermutedVector(mcVectors[currentMCVectorIndex], perm);
                    permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);
                    permutedRepresentativeIndices = icet::getPermutedVector(indicesOfRepresentativeSites, perm);

                    clusterVectorElement += _clusterSpace.evaluateClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first, permutedRepresentativeIndices) * elementsCountPair.second;
                }
            }

            /// This is the multiplicity one would have gotten during a full cluster vector calculation and is needed as normalizing factor
            double realMultiplicity = (double)_clusterSpace._sitePermutations[i][currentMCVectorIndex].size() * (double)_clusterSpace._orbitList._orbits[i].size() / (double)_clusterSpace._primitiveStructure.size();
            clusterVectorElement /= ((double)realMultiplicity * (double)_supercell.size());
            clusterVector.push_back(clusterVectorElement);
        }
    }
    return clusterVector;
}

/**
@details Calculate change in cluster vector upon change in occupation on one site
@param occupationsBefore the occupation vector for the supercell before the flip
@param flipIndex the index in the supercell where occupation has changed
@param newOccupation new atomic number on site index
*/
std::vector<double> ClusterExpansionCalculator::getClusterVectorChange(const std::vector<int> &occupationsBefore,
                                                                       int flipIndex,
                                                                       int newOccupation)
{
    _supercell.setAtomicNumbers(occupationsBefore);
    _clusterCounts.reset();
    if (occupationsBefore.size() != _supercell.size())
    {
        throw std::runtime_error("Input occupations and internal supercell structure mismatch in size (ClusterExpansionCalculator::getClusterVectorChange)");
    }
    if (flipIndex >= _supercell.size())
    {
        throw std::runtime_error("flipIndex larger than the length of the structure (ClusterExpansionCalculator::getClusterVectorChange)");
    }

    // do not sort the clusters
    bool keepOrder = true;

    // Count the clusters in the order as in the equivalent clusters
    // since these clusters are already in the permuted order
    bool permuteSites = false;

    // Get one of the translated orbitlists
    _translatedOrbitList = _localOrbitlists[_indexToOffset[flipIndex]];

    // Remove sites not containing the local index
    if (_clusterSpace._primitiveStructure.size() > 1)
    {
        // true meaning we only look at zero offset sites
        _translatedOrbitList.removeClustersWithoutIndex(flipIndex, true);
    }

    // Count clusters and get cluster count map
    _clusterCounts.countOrbitListChange(_supercell, flipIndex, newOccupation, _translatedOrbitList, keepOrder, permuteSites, -1, flipIndex);

    return _occupyClusterVector(0.0);
}

/**
@details This constructs a cluster vector that only includes clusters that contain the input index.
@param occupations the occupation vector for the supercell
@param index the local index of the supercell
*/
std::vector<double> ClusterExpansionCalculator::getLocalClusterVector(const std::vector<int> &occupations,
                                                                      int index)
{
    _supercell.setAtomicNumbers(occupations);
    _clusterCounts.reset();
    if (occupations.size() != _supercell.size())
    {
        throw std::runtime_error("Input occupations and internal supercell structure mismatch in size (ClusterExpansionCalculator::getLocalClusterVector)");
    }

    // do not sort the clusters
    bool keepOrder = true;

    // Count the clusters in the order as in the equivalent clusters
    // since these clusters are already in the permuted order
    bool permuteSites = false;

    // Get one of the translated orbitlists
    _translatedOrbitList = _localOrbitlists[_indexToOffset[index]];

    // Remove sites not containing the local index
    if (_clusterSpace._primitiveStructure.size() > 1)
    {
        // true meaning we only look at zero offset sites
        _translatedOrbitList.removeClustersWithoutIndex(index, true);
    }

    // Count clusters and get cluster count map
    _clusterCounts.countOrbitList(_supercell, _translatedOrbitList, keepOrder, permuteSites, -1, index);

    return _occupyClusterVector(1.0 / _supercell.size());
}

std::vector<double> ClusterExpansionCalculator::getClusterVector(const std::vector<int> &occupations)
{
    _supercell.setAtomicNumbers(occupations);
    _clusterCounts.reset();
    if (occupations.size() != _supercell.size())
    {
        throw std::runtime_error("Input occupations and internal supercell structure mismatch in size (ClusterExpansionCalculator::getLocalClusterVector)");
    }

    // do not sort the clusters
    bool keepOrder = true;

    /// Permute the sites
    bool permuteSites = true;

    // Count clusters and get cluster count map
    _clusterCounts.countOrbitList(_supercell, _fullOrbitList, keepOrder, permuteSites);

    return _occupyClusterVector(1.0);
}
