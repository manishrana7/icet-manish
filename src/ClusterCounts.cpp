#include "ClusterCounts.hpp"

/**
 @details Will count the vectors in latticeSites and assuming these sets of sites are represented by the cluster 'cluster'.
 @param structure the structure that will have its clusters counted
 @param latticeSites A group of sites, represented by 'cluster', that will be counted
 @param keepOrder if true the order of the sites will stay the same otherwise the vector of species being counted will be sorted
 @param siteIndexForDoubleCountCorrection In small supercells, clusters may include both a site and its periodic image.
                                      This argument can be used to avoid double counting in such cases; clusters
                                      in which a site with this index occurs more than once will only be counted
                                      with a factor 1 / n, where n is the number of occurences of this index.
                                      By default (siteIndexForDoubleCountCorrection = -1) no such correction is done.
*/
std::map<std::vector<int>, double> countClusters(const Structure &structure,
                                                 const std::vector<std::vector<LatticeSite>> &latticeSites,
                                                 bool keepOrder,
                                                 int siteIndexForDoubleCountCorrection)
{
    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elements(latticeSites[0].size());
    for (const auto &sites : latticeSites)
    {
        for (size_t i = 0; i < sites.size(); i++)
        {
            elements[i] = structure._atomicNumbers.at(sites[i].index());
        }
        if (!keepOrder)
        {
            std::sort(elements.begin(), elements.end());
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
        tmpCounts[elements] += unit;
    }
    return tmpCounts;
}

/**
 @details Count the change on occupation of the clusters of an orbit.
 @param structure the structure that will have its clusters counted
 @param flipIndex index of site that has been flipped
 @param newOccupation new atomic number of site that has been flipped
 @param latticeSites A group of sites, represented by 'cluster', that will be counted
 @param keepOrder if true the order of the sites will stay the same otherwise the vector of species being counted will be sorted
 @param siteIndexForDoubleCountCorrection In small supercells, clusters may include both a site and its periodic image.
                                      This argument can be used to avoid double counting in such cases; clusters
                                      in which a site with this index occurs more than once will only be counted
                                      with a factor 1 / n, where n is the number of occurences of this index.
                                      By default (siteIndexForDoubleCountCorrection = -1) no such correction is done.
*/
std::map<std::vector<int>, double> countClusterChanges(const Structure &structure,
                                                       const int flipIndex,
                                                       const int newOccupation,
                                                       const std::vector<std::vector<LatticeSite>> &latticeSites,
                                                       bool keepOrder,
                                                       int siteIndexForDoubleCountCorrection)
{
    std::map<std::vector<int>, double> tmpCounts;
    std::vector<int> elementsOld(latticeSites[0].size());
    std::vector<int> elementsNew(latticeSites[0].size());
    int siteIndex;
    int occupation;
    for (const auto &sites : latticeSites)
    {
        for (size_t i = 0; i < sites.size(); i++)
        {
            siteIndex = sites[i].index();
            occupation = structure._atomicNumbers.at(siteIndex);
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
        if (!keepOrder)
        {
            std::sort(elementsOld.begin(), elementsOld.end());
            std::sort(elementsNew.begin(), elementsNew.end());
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

/**
 @brief Counts the clusters in the input structure.
 @param structure input configuration
 @param orbitList orbit list
 @param keepOrder if true do not reorder clusters before comparison (i.e., ABC != ACB)
 @param permuteSites if true the sites will be permuted according to the correspondin permutations in the orbit
 @param maxOrbit include only orbits with indices smaller than this (by default all orbits are included)
 @param siteIndexForDoubleCountCorrection In small supercells, clusters may include both a site and its periodic image.
                                      This argument can be used to avoid double counting in such cases; clusters
                                      in which a site with this index occurs more than once will only be counted
                                      with a factor 1 / n, where n is the number of occurences of this index.
                                      By default (siteIndexForDoubleCountCorrection = -1) no such correction is done.
*/
std::unordered_map<Cluster, std::map<std::vector<int>, double>> countOrbitList(const Structure &structure, const OrbitList &orbitList,
                                                                               bool keepOrder, bool permuteSites, int maxOrbit,
                                                                               int siteIndexForDoubleCountCorrection)
{
    std::unordered_map<Cluster, std::map<std::vector<int>, double>> clusterCounts;
    if (maxOrbit == -1)
    {
        maxOrbit = orbitList.size();
    }
    for (size_t i = 0; i < maxOrbit; i++)
    {
        Cluster representativeCluster = orbitList._orbits[i].getRepresentativeCluster();
        representativeCluster.setTag(i);
        std::map<std::vector<int>, double> partialCounts;
        if (permuteSites && keepOrder && representativeCluster.order() != 1)
        {
            partialCounts = countClusters(structure, orbitList.getOrbit(i).getPermutedEquivalentClusters(), keepOrder, siteIndexForDoubleCountCorrection);
        }
        else
        {
            partialCounts = countClusters(structure, orbitList._orbits[i]._equivalentClusters, keepOrder, siteIndexForDoubleCountCorrection);
        }
        // Now add counts to the "master" _clusterCounts
        for (auto count : partialCounts)
        {
            clusterCounts[representativeCluster][count.first] += count.second;
        }
    }
    return clusterCounts;
}

/**
 @brief Counts the clusters in the input structure.
 @param structure input configuration
 @param flipIndex index of site that has been flipped
 @param newOccupation new atomic number of site that has been flipped
 @param orbitList orbit list
 @param keepOrder if true do not reorder clusters before comparison (i.e., ABC != ACB)
 @param maxOrbit include only orbits with indices smaller than this (by default all orbits are included)
 @param siteIndexForDoubleCountCorrection In small supercells, clusters may include both a site and its periodic image.
                                      This argument can be used to avoid double counting in such cases; clusters
                                      in which a site with this index occurs more than once will only be counted
                                      with a factor 1 / n, where n is the number of occurences of this index.
                                      By default (siteIndexForDoubleCountCorrection = -1) no such correction is done.
*/
std::unordered_map<Cluster, std::map<std::vector<int>, double>> countOrbitListChange(const Structure &structure,
                                                                                     const int flipIndex,
                                                                                     const int newOccupation,
                                                                                     const OrbitList &orbitList,
                                                                                     bool keepOrder,
                                                                                     int maxOrbit,
                                                                                     int siteIndexForDoubleCountCorrection)
{
    std::unordered_map<Cluster, std::map<std::vector<int>, double>> clusterCounts;
    if (maxOrbit == -1)
    {
        maxOrbit = orbitList.size();
    }
    for (size_t i = 0; i < maxOrbit; i++)
    {
        Cluster representativeCluster = orbitList._orbits[i].getRepresentativeCluster();
        representativeCluster.setTag(i);
        // Here we rely on _equivalentClusters being pre-permuted, as is done in the ClusterExpansionCalculator
        std::map<std::vector<int>, double> partialCounts = countClusterChanges(structure, flipIndex, newOccupation, orbitList._orbits[i]._equivalentClusters, keepOrder, siteIndexForDoubleCountCorrection);

        // Now add counts to the "master" _clusterCounts
        for (auto count : partialCounts)
        {
            clusterCounts[representativeCluster][count.first] += count.second;
        }
    }
    return clusterCounts;
}
