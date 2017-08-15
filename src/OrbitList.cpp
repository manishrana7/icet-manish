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

        }
    }

}