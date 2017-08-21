#include "OrbitList.hpp"

/**
@TODO: Think about adding a string tag here to keep track of different orbitlists
*/
OrbitList::OrbitList()
{
    //Empty constructor
}

///Construct orbitlist from mbnl and structure
OrbitList::OrbitList(const ManybodyNeighborlist &mbnl, const Structure &structure)
{
    std::unordered_map<Cluster, int> clusterIndexMap;
    for (size_t i = 0; i < mbnl.getNumberOfSites(); i++)
    {
        //special case for singlet
        if (mbnl.getNumberOfSites(i) == 0)
        {
            std::vector<LatticeNeighbor> sites = mbnl.getSites(i, 0);
            Cluster cluster = Cluster(structure, sites);
            addClusterToOrbitlist(cluster, sites,clusterIndexMap);
        }

        for (size_t j = 0; j < mbnl.getNumberOfSites(i); j++)
        {
            std::vector<LatticeNeighbor> sites = mbnl.getSites(i, j);
            Cluster cluster = Cluster(structure, sites);
            addClusterToOrbitlist(cluster, sites,clusterIndexMap);
        }
    }
}


///add cluster to orbitlist, if cluster exists add sites if not create a new orbit
void OrbitList::addClusterToOrbitlist(const Cluster &cluster, const std::vector<LatticeNeighbor> &sites, std::unordered_map<Cluster, int> &clusterIndexMap)
{
    int orbitNumber = findOrbit(cluster, clusterIndexMap);
    if (orbitNumber == -1)
    {
        Orbit newOrbit = Orbit(cluster);
        addOrbit(newOrbit);
        //add to back ( assuming addOrbit does not sort orbitlist )
        _orbitList.back().addEquivalentSites(sites);
        clusterIndexMap[cluster] = _orbitList.size() - 1;
    }
    else
    {
        _orbitList[orbitNumber].addEquivalentSites(sites);
    }
}

/**
Returns the orbit for which "cluster" is the representative cluster

returns -1 if it nothing is found
*/
int OrbitList::findOrbit(const Cluster &cluster) const
{
    for (size_t i = 0; i < _orbitList.size(); i++)
    {
        if (_orbitList[i].getRepresentativeCluster() == cluster)
        {
            return i;
        }
    }
    return -1;
}

/**
Returns the orbit for which "cluster" is the representative cluster

returns -1 if it nothing is found
*/
int OrbitList::findOrbit(const Cluster &cluster, const std::unordered_map<Cluster, int> &clusterIndexMap) const
{
    auto search = clusterIndexMap.find(cluster);
    if (search != clusterIndexMap.end())
    {
        return search->second;
    }
    else
    {
        return -1;
    }
}
