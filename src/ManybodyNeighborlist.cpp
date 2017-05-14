#include "ManybodyNeighborlist.hpp"

std::vector<std::pair<int, Vector3d>> ManybodyNeighborlist::build(const Neighborlist &nl, int index, int order)
{
    std::vector<std::pair<int, Vector3d>> manybodyNeighborIndex;
    auto Ni = nl.getNeighbors(index);
    int numberOfSites = nl.size();

    std::vector<int> neighborIndices;
    neighborIndices.push_back(index);

    for (int c = 2; c < order; c++)
    {
        goDeeper(nl,neighborIndices,manybodyNeighborIndex, Ni,index);       
    }
}


void ManybodyNeighborlist::goDeeper(const Neighborlist &nl, std::vector<int> neighborIndices,
                                 std::vector<std::pair<int, Vector3d>> manybodyNeighborIndex, std::vector<std::pair<int, Vector3d>> Ni, int index )
{

    for (int i = 0; i < nl.size(); i++)
        {
            if (nl.isNeighbor(index, i))
            {
                neighborIndices.push_back(i);
                auto intersection = getIntersection(Ni, nl.getNeighbors(i));
            }
        }
}