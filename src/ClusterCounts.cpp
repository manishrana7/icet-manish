#include "ClusterCounts.hpp"

/**
@details Will count the vectors in latticeSites and assuming these sets of sites are represented by the cluster 'cluster'.
@param structure the structure that will have its clusters counted
@param latticeSites A group of sites, represented by 'cluster', that will be counted
@param cluster A cluster used as identification on what sites the clusters belong to
@param keepOrder if true the order of the sites will stay the same otherwise the vector of species being counted will be sorted
*/
void ClusterCounts::count(const Structure &structure, const std::vector<std::vector<LatticeSite>> &latticeSites,
                          const Cluster &cluster, bool keepOrder)
{
    std::vector<int> elements(latticeSites[0].size());
    for (const auto &sites : latticeSites)
    {
        for (size_t i = 0; i < sites.size(); i++)
        {
            elements[i] = structure._atomicNumbers.at(sites[i].index());
        }
        countCluster(cluster, elements, keepOrder);
    }
}

/// Counts cluster only through this function.
void ClusterCounts::countCluster(const Cluster &cluster, const std::vector<int> &elements, bool keepOrder, int unit)
{
    if (keepOrder)
    {
        _clusterCounts[cluster][elements] += unit;
    }
    else
    {
        std::vector<int> sortedElements = elements;
        std::sort(sortedElements.begin(), sortedElements.end());
        _clusterCounts[cluster][sortedElements] += unit;
    }
}

/**
 @brief Counts the clusters in the input structure.
 @param structure input configuration
 @param orbitList orbit list
 @param keepOrder if true do not reorder clusters before comparison (i.e., ABC != ACB)
 @param permuteSites if true the sites will be permuted according to the correspondin permutations in the orbit
 @param maxOrbit include only orbits with indices smaller than this (by default all orbits are included)
*/
void ClusterCounts::countOrbitList(const Structure &structure, const OrbitList &orbitList, bool keepOrder, bool permuteSites, int maxOrbit)
{
    if (maxOrbit == -1)
    {
        maxOrbit = orbitList.size();
    }
    for (size_t i = 0; i < maxOrbit; i++)
    {
        Cluster representativeCluster = orbitList._orbits[i].getRepresentativeCluster();
        representativeCluster.setTag(i);
        if (permuteSites && keepOrder && representativeCluster.order() != 1)
        {
            count(structure, orbitList.getOrbit(i).getPermutedEquivalentClusters(), representativeCluster, keepOrder);
        }
        else if (!permuteSites && keepOrder && representativeCluster.order() != 1)
        {
            count(structure, orbitList._orbits[i]._equivalentClusters, representativeCluster, keepOrder);
        }
        else
        {
            count(structure, orbitList._orbits[i]._equivalentClusters, representativeCluster, keepOrder);
        }
    }
}

/**
@details Will count the vectors in latticeSites and assuming these sets of sites are represented by the cluster 'cluster'.
@param structure the structure that will have its clusters counted
@param flipIndex index of site that has been flipped
@param newOccupation new atomic number of site that has been flipped
@param latticeSites A group of sites, represented by 'cluster', that will be counted
@param cluster A cluster used as identification on what sites the clusters belong to
@param keepOrder if true the order of the sites will stay the same otherwise the vector of species being counted will be sorted
*/
void ClusterCounts::countChange(const Structure &structure,
                                const int flipIndex,
                                const int newOccupation,
                                const std::vector<std::vector<LatticeSite>> &latticeSites,
                                const Cluster &cluster, bool keepOrder)
{
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
        // Now count the change in elements, -1 for every combination of elements
        // that disappear with the new occupation, +1 for every combination that
        // appeared
        countCluster(cluster, elementsOld, keepOrder, -1);
        countCluster(cluster, elementsNew, keepOrder, 1);
    }
}

/**
 @brief Counts the clusters in the input structure.
 @param structure input configuration
 @param flipIndex index of site that has been flipped
 @param newOccupation new atomic number of site that has been flipped
 @param orbitList orbit list
 @param keepOrder if true do not reorder clusters before comparison (i.e., ABC != ACB)
 @param permuteSites if true the sites will be permuted according to the correspondin permutations in the orbit
 @param maxOrbit include only orbits with indices smaller than this (by default all orbits are included)
*/
void ClusterCounts::countOrbitListChange(const Structure &structure,
                                         const int flipIndex,
                                         const int newOccupation,
                                         const OrbitList &orbitList,
                                         bool keepOrder, bool permuteSites, int maxOrbit)
{
    if (maxOrbit == -1)
    {
        maxOrbit = orbitList.size();
    }
    for (size_t i = 0; i < maxOrbit; i++)
    {
        Cluster representativeCluster = orbitList._orbits[i].getRepresentativeCluster();
        representativeCluster.setTag(i);
        if (permuteSites && keepOrder && representativeCluster.order() != 1)
        {
            countChange(structure, flipIndex, newOccupation, orbitList.getOrbit(i).getPermutedEquivalentClusters(), representativeCluster, keepOrder);
        }
        else
        {
            countChange(structure, flipIndex, newOccupation, orbitList._orbits[i]._equivalentClusters, representativeCluster, keepOrder);
        }
    }
}
