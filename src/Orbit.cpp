#include "Orbit.hpp"

/**
@details This constructor creates an orbit from a set of equivalent clusters and a structure.
         Note that the sites of each cluster need to be ordered in a consistent manner
         (This ordering is enforced when Orbit objects are created via an OrbitList object)
@param clusters
    The clusters that together make up the orbit. Normally, these should be symmetry equivalent,
    but note that this needs to be ensured before constructing the orbit; it is not enforced by
    the constructor.
@param allowedClusterPermutations
    Allowed permutations for this orbit; e.g., if 0, 2, 1 is in this set
    then 0, 1, 0 is the same multi-component vector as 0, 0, 1
**/
Orbit::Orbit(const std::vector<Cluster> clusters,
             const std::set<std::vector<int>> allowedClusterPermutations)
{
    _representativeCluster = clusters[0];
    _clusters = clusters;
    _allowedClusterPermutations = allowedClusterPermutations;

    // Sort clusters according the coordinates of their sites.
    // This is done only to achieve a reproducible ordering, which under some
    // circumstances matters when sorting the orbit list (and thereby the
    // cluster vector).
    std::sort(_clusters.begin(), _clusters.end());
}

/**
@brief Add a cluster to the orbit.
@details Note that this function only appends the new cluster to the end without resorting.
@param latticeSiteGroup Cluster to be added represented by a group of lattice sites
*/
void Orbit::addCluster(const Cluster &cluster)
{
    _clusters.push_back(cluster);
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
@details Check if this orbit contains a cluster in its list of clusters.
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

    for (size_t i = 0; i < _clusters.size(); i++)
    {
        auto sites = _clusters[i].getLatticeSites();

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
std::map<std::vector<int>, double> Orbit::countClusters(const std::shared_ptr<Structure> structure,
                                                        int siteIndexForDoubleCountingCorrection) const
{
    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elements(order());
    for (const auto &cluster : _clusters)
    {

        // If we apply the double counting correction for some site we need
        // to ensure that we only count clusters that include this site.
        if (siteIndexForDoubleCountingCorrection < 0 || cluster.isSiteIndexIncludedWithZeroOffset(siteIndexForDoubleCountingCorrection))
        {
            const std::vector<LatticeSite> &sites = cluster.getLatticeSites();
            for (size_t i = 0; i < sites.size(); i++)
            {
                elements[i] = structure->getAtomicNumbers().at(sites[i].index());
            }
            // If the current atom (siteIndexForDoubleCountingCorrection) occurs more than once,
            // we risk double counting it if we calculate a local cluster vector. To avoid this,
            // we count the clusters in units of 1 / n, where n is the number of occurences of
            // the present atom in the cluster.
            double unit = 1;
            if (siteIndexForDoubleCountingCorrection > -1)
            {
                unit /= cluster.countOccurencesOfSiteIndex(siteIndexForDoubleCountingCorrection);
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
std::map<std::vector<int>, double> Orbit::countClusterChanges(const std::shared_ptr<Structure> structure,
                                                              const int flipIndex,
                                                              const int newOccupation) const
{
    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elementsOld(order());
    std::vector<int> elementsNew(order());
    int siteIndex;
    int occupation;

    for (const auto &cluster : _clusters)
    {
        // Only count clusters where site flipIndex is included (with zero offset)
        if (cluster.isSiteIndexIncludedWithZeroOffset(flipIndex))
        {
            const std::vector<LatticeSite> &sites = cluster.getLatticeSites();
            for (size_t i = 0; i < sites.size(); i++)
            {
                siteIndex = sites[i].index();
                occupation = structure->getAtomicNumbers().at(siteIndex);
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
            // If the current site (flipIndex) occurs more than once,
            // we risk double counting it. To avoid this, we count the clusters in units of
            // 1 / n, where n is the number of occurrences of the present atom in the cluster.
            double unit = 1.0 / (double)cluster.countOccurencesOfSiteIndex(flipIndex);

            // The old cluster has disappeared and we got elementNew instead
            tmpCounts[elementsOld] -= unit;
            tmpCounts[elementsNew] += unit;
        }
    }
    return tmpCounts;
}

/**
**/
void Orbit::translate(const Vector3d &cellOffset) 
{
    _representativeCluster.translate(cellOffset);
    for (auto &cluster : _clusters) {
        cluster.translate(cellOffset);
    }
}

/**
@brief Transforms the clusters to a new cell
@details
    Each cluster in the orbit consists of a vector of lattice sites, and
    these sites are defined in relation to a specific atomic structure
    (typically a primitive structure), and a pointer to this structure
    is stored in each cluster. This function redefines the sites
    such that they refer to a new cell (typically a supercell).
@param supercell The new atomic structure
@param cellOffset
    Offset to be applied to sites before transformation to supercell.
    This offset is specified in terms of the old (primitive) structure. 
@param primToSuperMap
    Map from lattice site referring to old structure to lattice site
    referring to the new structure. This map will successivelly be
    populated when executing the function, and is only used for
    reasons of performance.
@fractionalPositionTolerance 
**/
void Orbit::transformToSupercell(std::shared_ptr<Structure> supercellPtr,
                                 std::unordered_map<LatticeSite, LatticeSite> &primToSuperMap,
                                 const double fractionalPositionTolerance)
{   
    _representativeCluster.transformToSupercell(supercellPtr, primToSuperMap, fractionalPositionTolerance);
    for (auto &cluster : _clusters)
    {
        cluster.transformToSupercell(supercellPtr, primToSuperMap, fractionalPositionTolerance);
    }
}

namespace std
{

    /// Stream operator.
    ostream &operator<<(ostream &os, const Orbit &orbit)
    {
        for (const auto cluster : orbit.getClusters())
        {
            os << "  ";
            for (const auto site : cluster.getLatticeSites())
            {
                os << " " << site;
            }
        }
        return os;
    }

}
