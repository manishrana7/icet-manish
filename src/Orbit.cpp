#include "Orbit.hpp"

/**
@details This constructor creates an orbit from a set of equivalent clusters and a structure.
         Note that the sites of each cluster need to be ordered in a consistent manner
         (this ordering is enforced when Orbits are created from the OrbitList class)
@param equivalentClusters The clusters that together make up the orbit.
@param structure structure that this orbit relates to
@param allowedClusterPermutations
    allowed permutations for this orbit; e.g., if 0,2,1 is in this set
    then 0,1,0 is the same MC vector as 0,0,1
**/
Orbit::Orbit(const Structure &structure,
             const std::vector<std::vector<LatticeSite>> equivalentClusters,
             const std::set<std::vector<int>> allowedClusterPermutations)
{
    _representativeCluster = Cluster(structure, equivalentClusters[0]);
    _equivalentClusters = equivalentClusters;
    _allowedClusterPermutations = allowedClusterPermutations;
    sort();
}

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
@brief Check whether a site is included in a cluster
@details A cluster will count as included if index is among the lattice sites
         that have a zero offset.
@param index Index of site to check whether it is included
@param cluster
    Vector of LatticeSite, at least one of which must contain index for this
    function to return true
*/
bool Orbit::isSiteIncluded(int index, const std::vector<LatticeSite> &cluster) const
{
    return std::any_of(cluster.begin(), cluster.end(), [=](const LatticeSite &ls)
                       { return ls.index() == index && ls.unitcellOffset().norm() < 1e-4; });
}

/**
 @brief Count the occupations of the clusters in this orbit.
 @details
    Note that the orderings of the sites in the clusters matter, meaning,
    for example, that (47, 79) will be counted separately from (79, 47)
    (here 47 and 79 are atomic numbers).
 @param structure the structure that will have its clusters counted
 @param siteIndexForDoubleCountingCorrection
   In small supercells, clusters may include both a site and its periodic image.
   In such cases this argument can be used to avoid double counting.
   Clusters in which a site with this index occurs more than once will only be counted with
   a factor 1/n, where n is the number of occurrences of this index. By default
   (i.e. siteIndexForDoubleCountingCorrection = -1) no such correction is applied.
*/
std::map<std::vector<int>, double> Orbit::countClusters(const Structure &structure,
                                                        int siteIndexForDoubleCountingCorrection) const
{
    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elements(order());
    for (const auto &sites : _equivalentClusters)
    {
        // If we apply the double counting correction for some site we need
        // to ensure that we only count clusters that include this site.
        if (siteIndexForDoubleCountingCorrection < 0 || isSiteIncluded(siteIndexForDoubleCountingCorrection, sites))
        {
            for (size_t i = 0; i < sites.size(); i++)
            {
                elements[i] = structure.getAtomicNumbers().at(sites[i].index());
            }
            double unit = 1;
            // If the current atom (siteIndexForDoubleCountingCorrection) occurs more than once,
            // we risk double counting it if we calculate a change in cluster vector or
            // a local cluster vector. To avoid this, we count the clusters in units of
            // 1 / n, where n is the number of occurences of the present atom in the cluster.
            if (siteIndexForDoubleCountingCorrection > -1)
            {
                unit /= (double)std::count_if(sites.begin(), sites.end(), [=](LatticeSite ls)
                                              { return ls.index() == siteIndexForDoubleCountingCorrection; });
            }
            tmpCounts[elements] += unit;
        }
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
*/
std::map<std::vector<int>, double> Orbit::countClusterChanges(const Structure &structure,
                                                              const int flipIndex,
                                                              const int newOccupation) const
{
    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elementsOld(order());
    std::vector<int> elementsNew(order());
    int siteIndex;
    int occupation;

    for (const auto &sites : _equivalentClusters)
    {
        // Only count clusters where site flipIndex is included (with zero offset)
        if (isSiteIncluded(flipIndex, sites))
        {
            for (size_t i = 0; i < sites.size(); i++)
            {
                siteIndex = sites[i].index();
                occupation = structure.getAtomicNumbers().at(siteIndex);
                elementsOld[i] = occupation;

                // If the present site is the one that was changed,
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
            // If the current site (flipIndex) occurs more than once,
            // we risk double counting it. To avoid this, we count the clusters in units of
            // 1 / n, where n is the number of occurrences of the present atom in the cluster.
            unit /= (double)std::count_if(sites.begin(), sites.end(), [=](LatticeSite ls)
                                          { return ls.index() == flipIndex; });
            // The old cluster has disappeared and we got elementNew instead
            tmpCounts[elementsOld] -= unit;
            tmpCounts[elementsNew] += unit;
        }
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
