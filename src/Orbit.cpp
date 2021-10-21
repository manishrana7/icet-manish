#include "Orbit.hpp"

/**
@param latticeSiteGroup cluster to be added represented by a group of lattice site
@param sort_orbit if True the orbit will be sorted
*/
void Orbit::addEquivalentCluster(const std::vector<LatticeSite> &latticeSiteGroup, bool sort_orbit)
{
    _equivalentClusters.push_back(latticeSiteGroup);
    if (sort_orbit)
    {
        sort();
    }
}

/**
@param latticeSiteGroups list of cluster to be added, each represented by a group of lattice sites
@param sort_orbit if True the orbit will be sorted
*/
void Orbit::addEquivalentClusters(const std::vector<std::vector<LatticeSite>> &latticeSiteGroups, bool sort_orbit)
{
    _equivalentClusters.insert(_equivalentClusters.end(), latticeSiteGroups.begin(), latticeSiteGroups.end());
    if (sort_orbit)
    {
        sort();
    }
}

/**
@param index cluster index
*/
std::vector<LatticeSite> Orbit::getClusterByIndex(unsigned int index) const
{
    if (index >= _equivalentClusters.size())
    {
        throw std::out_of_range("Index out of range (Orbit::getClusterByIndex)");
    }
    return _equivalentClusters[index];
}

/**
@param index cluster index
*/
std::vector<LatticeSite> Orbit::getPermutedClusterByIndex(unsigned int index) const
{
    if (index >= _equivalentClusters.size())
    {
        std::ostringstream msg;
        msg << "Index out of range (Orbit::getPermutedClusterByIndex).\n";
        msg << " index: " << index << "\n";
        msg << " size(_equivalentClusters): " << _equivalentClusters.size() << "\n";
        throw std::out_of_range(msg.str());
    }
    if (index >= _equivalentClusterPermutations.size())
    {
        std::ostringstream msg;
        msg << "Index out of range (Orbit::getPermutedClusterByIndex).\n";
        msg << " index: " << index << "\n";
        msg << " size(_equivalentClusterPermutations): " << _equivalentClusterPermutations.size();
        throw std::out_of_range(msg.str());
    }
    return icet::getPermutedVector<LatticeSite>(_equivalentClusters[index], _equivalentClusterPermutations[index]);
}

/**
@brief Returns all permuted equivalent sites.
*/
std::vector<std::vector<LatticeSite>> Orbit::getPermutedEquivalentClusters() const
{
    std::vector<std::vector<LatticeSite>> permutedSites(_equivalentClusters.size());
    for (size_t i = 0; i < _equivalentClusters.size(); i++)
    {
        permutedSites[i] = getPermutedClusterByIndex(i);
    }
    return permutedSites;
}

/**
  @param Mi_local list of the number of allowed species per site
 */
std::vector<std::vector<int>> Orbit::getMultiComponentVectors(const std::vector<int> &Mi_local) const
{
    if (std::any_of(Mi_local.begin(), Mi_local.end(), [](const int i)
                    { return i < 2; }))
    {
        std::vector<std::vector<int>> emptyVector;
        return emptyVector;
    }
    auto allMCVectors = getAllPossibleMultiComponentVectorPermutations(Mi_local);
    std::sort(allMCVectors.begin(), allMCVectors.end());
    std::vector<std::vector<int>> distinctMCVectors;
    for (const auto &mcVector : allMCVectors)
    {
        std::vector<std::vector<int>> permutedMCVectors;
        for (const auto &allowedPermutation : _allowedClusterPermutations)
        {
            permutedMCVectors.push_back(icet::getPermutedVector<int>(mcVector, allowedPermutation));
        }
        std::sort(permutedMCVectors.begin(), permutedMCVectors.end());

        // if not any of the vectors in permutedMCVectors exist in distinctMCVectors
        if (!std::any_of(permutedMCVectors.begin(), permutedMCVectors.end(), [&](const std::vector<int> &permMcVector)
                         { return !(std::find(distinctMCVectors.begin(), distinctMCVectors.end(), permMcVector) == distinctMCVectors.end()); }))
        {
            distinctMCVectors.push_back(mcVector);
        }
    }
    return distinctMCVectors;
}

/// Similar to get all permutations but needs to be filtered through the number of allowed species
std::vector<std::vector<int>> Orbit::getAllPossibleMultiComponentVectorPermutations(const std::vector<int> &Mi_local) const
{

    std::vector<std::vector<int>> cartesianFactors(Mi_local.size());
    for (size_t i = 0; i < Mi_local.size(); i++)
    {
        for (int j = 0; j < Mi_local[i] - 1; j++) // N.B minus one so a binary only get one cluster function
        {
            cartesianFactors[i].push_back(j);
        }
    }

    std::vector<std::vector<int>> allPossibleMCPermutations;
    std::vector<int> firstVector(Mi_local.size(), 0);

    do
    {
        allPossibleMCPermutations.push_back(firstVector);
    } while (icet::nextCartesianProduct(cartesianFactors, firstVector));

    return allPossibleMCPermutations;
}

/**
@details Check if this orbit contains a cluster in its list of equivalent clusters.
@param cluster cluster (represented by a list of sites) to look for
@param sorted if true the order of sites in the cluster is irrelevant
@returns true if the cluster is present in the orbit
**/
bool Orbit::contains(const std::vector<LatticeSite> cluster, bool sorted) const
{
    auto clusterCopy = cluster;
    if (sorted)
    {
        std::sort(clusterCopy.begin(), clusterCopy.end());
    }

    for (size_t i = 0; i < _equivalentClusters.size(); i++)
    {
        auto sites = _equivalentClusters[i];

        // compare the sorted sites
        if (sorted)
        {
            std::sort(sites.begin(), sites.end());
        }

        if (sites == clusterCopy)
        {
            return true;
        }
    }
    return false;
}

/**
@param index site index
@param onlyConsiderSitesWithZeroOffset if true only consider sites with zero offset
*/
void Orbit::removeClustersWithoutIndex(const size_t index, bool onlyConsiderSitesWithZeroOffset)
{
    for (int i = _equivalentClusters.size() - 1; i >= 0; i--)
    {
        if (onlyConsiderSitesWithZeroOffset)
        {
            if (std::none_of(_equivalentClusters[i].begin(), _equivalentClusters[i].end(), [=](LatticeSite &ls)
                             { return ls.index() == index && ls.unitcellOffset().norm() < 1e-4; }))
            {
                _equivalentClusters.erase(_equivalentClusters.begin() + i);
                _equivalentClusterPermutations.erase(_equivalentClusterPermutations.begin() + i);
            }
        }
        else
        {
            if (std::none_of(_equivalentClusters[i].begin(), _equivalentClusters[i].end(), [=](LatticeSite &ls)
                             { return ls.index() == index; }))
            {
                _equivalentClusters.erase(_equivalentClusters.begin() + i);
                _equivalentClusterPermutations.erase(_equivalentClusterPermutations.begin() + i);
            }
        }
    }
}

/**
@param cluster cluster to be removed (represented by a list of sites); the order of sites is irrelevant
*/
void Orbit::removeCluster(std::vector<LatticeSite> cluster)
{

    std::sort(cluster.begin(), cluster.end());
    for (size_t i = 0; i < _equivalentClusters.size(); i++)
    {
        auto sites = _equivalentClusters[i];

        // compare the sorted sites
        std::sort(sites.begin(), sites.end());

        if (sites == cluster)
        {
            _equivalentClusters.erase(_equivalentClusters.begin() + i);
            _equivalentClusterPermutations.erase(_equivalentClusterPermutations.begin() + i);
            return;
        }
    }
    throw std::runtime_error("Did not find any matching clusters (Orbit::removeCluster)");
}

/**
 @brief Count the occupations of the clusters in this orbit.
 @details
    Note that the orderings of the sites in the clusters matter, meaning,
    for example, that (47, 79) will be counted separately from (79, 47)
    (here 47 and 79 are atomic numbers).
 @param structure the structure that will have its clusters counted
 @param siteIndexForDoubleCountCorrection
   In small supercells, clusters may include both a site and its periodic image.
   In such cases this argument can be used to avoid double counting.
   Clusters in which a site with this index occurs more than once will only be counted with
   a factor 1/n, where n is the number of occurrences of this index. By default
   (i.e. siteIndexForDoubleCountCorrection = -1) no such correction is applied.
 @param permuteClusters If true, permute clusters equivalent clusters before counting
*/
std::map<std::vector<int>, double> Orbit::countClusters(const Structure &structure,
                                                        int siteIndexForDoubleCountCorrection,
                                                        bool permuteClusters) const
{
    if (permuteClusters)
    {
        // In this case we could just loop over getPermutedEquivalentClusters() instead of
        // _equivalentClusters, but we then need to be very careful with performance.
        // Since this should never happen unless someone does something very specific,
        // we disallow it for now.
        throw std::runtime_error("countClusterChanges does not support counting of clusters that are not permuted (Orbit::coundClusterChanges)");
    }

    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elements(order());
    for (const auto &sites : _equivalentClusters)
    {
        for (size_t i = 0; i < sites.size(); i++)
        {
            elements[i] = structure.getAtomicNumbers().at(sites[i].index());
        }
        double unit = 1;
        // If the current atom (siteIndexForDoubleCountCorrection) occurs more than once,
        // we risk double counting it if we calculate a change in cluster vector or
        // a local cluster vector. To avoid this, we count the clusters in units of
        // 1 / n, where n is the number of occurences of the present atom in the cluster.
        if (siteIndexForDoubleCountCorrection > -1)
        {
            unit /= (double)std::count_if(sites.begin(), sites.end(), [=](LatticeSite ls)
                                          { return ls.index() == siteIndexForDoubleCountCorrection; });
        }
        tmpCounts[elements] += unit;
    }
    return tmpCounts;
}

/**
 @brief 
    Count the change in occupations of the clusters in this orbit caused
    by changing the chemical identity of one site.
 @details
    `structure` should contain the original occupations, and the change
    is defined by `flipIndex` (index of the site whose occupation
    changes) and `newOccupation` (the new atomic number on that site).
    Note that the orderings of the sites in the clusters matter, meaning,
    for example, that (47, 79) will be counted separately from (79, 47)
    (here 47 and 79 are atomic numbers).
 @param structure the structure for which to count clusters, with occupations before change
 @param flipIndex index of site that has been flipped
 @param newOccupation new atomic number of site that has been flipped
 @param siteIndexForDoubleCountCorrection
   In small supercells, clusters may include both a site and its periodic image.
   In such cases this argument can be used to avoid double counting.
   Clusters in which a site with this index occurs more than once will only be counted with
   a factor 1/n, where n is the number of occurrences of this index. By default
   (i.e. siteIndexForDoubleCountCorrection = -1) no such correction is applied.
 @param permuteClusters If true, permute clusters equivalent clusters before counting (not yet implemented, should normally never be done)
*/
std::map<std::vector<int>, double> Orbit::countClusterChanges(const Structure &structure,
                                                              const int flipIndex,
                                                              const int newOccupation,
                                                              int siteIndexForDoubleCountCorrection,
                                                              const bool permuteClusters) const
{
    if (permuteClusters)
    {
        // In this case we could just loop over getPermutedEquivalentClusters() instead of
        // _equivalentClusters, but we then need to be very careful with performance.
        // Since this should never happen unless someone does something very specific,
        // we disallow it for now.
        throw std::runtime_error("countClusterChanges does not support counting of clusters that are not permuted (Orbit::coundClusterChanges)");
    }

    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elementsOld(order());
    std::vector<int> elementsNew(order());
    int siteIndex;
    int occupation;

    for (const auto &sites : _equivalentClusters)
    {
        for (size_t i = 0; i < sites.size(); i++)
        {
            siteIndex = sites[i].index();
            occupation = structure.getAtomicNumbers().at(siteIndex);
            elementsOld[i] = occupation;

            // If the present site index is the one that was changed,
            // we need to use a different atomic number
            if (siteIndex == flipIndex)
            {
                elementsNew[i] = newOccupation;
            }
            else
            {
                elementsNew[i] = occupation;
            }
        }
        double unit = 1;
        // If the present atom (siteIndexForDoubleCountCorrection) occurs more than once,
        // we risk double counting it if we calculate a change in cluster vector or
        // a local cluster vector. To avoid this, we count the clusters in units of
        // 1 / n, where n is the number of occurences of the present atom in the cluster.
        if (siteIndexForDoubleCountCorrection > -1)
        {
            unit /= (double)std::count_if(sites.begin(), sites.end(), [=](LatticeSite ls)
                                          { return ls.index() == siteIndexForDoubleCountCorrection; });
        }
        // The old cluster has disappeared and we have gotten elementNew instead; count that
        tmpCounts[elementsOld] -= unit;
        tmpCounts[elementsNew] += unit;
    }
    return tmpCounts;
}

namespace std
{

    /// Stream operator.
    ostream &operator<<(ostream &os, const Orbit &orbit)
    {
        for (const auto cluster : orbit.getEquivalentClusters())
        {
            os << "  ";
            for (const auto site : cluster)
            {
                os << " " << site;
            }
        }
        return os;
    }

}
