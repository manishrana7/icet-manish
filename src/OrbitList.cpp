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
    for(size_t i=0; i< mbnl.getNumberOfSites(); i++)
    {
        for(size_t j=0; j<mbnl.getNumberOfSites(i); j++)
        {
            std::vector<LatticeNeighbor> sites = mbnl.getSites(i,j);
            Cluster cluster = Cluster(structure, sites);
            int orbitNumber = findOrbit(cluster);
            if(orbitNumber == -1)
            {
                Orbit newOrbit = Orbit(cluster);
                addOrbit(newOrbit);
            }
            else
            {
                _orbitList[orbitNumber].addEquivalentSites(sites);
            }
        }
    }

}

/**
Returns the orbit for which "cluster" is the representative cluster

returns -1 if it nothing is found
*/
int OrbitList::findOrbit(const Cluster &cluster) const
{
    for(size_t i = 0; i < _orbitList.size(); i++)
    {
        if(_orbitList[i].getRepresentativeCluster() == cluster)
        {
            return i;
        }
    }

    return -1;

}