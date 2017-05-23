#include "ManybodyNeighborlist.hpp"

/**
    This will use the neighborlist to combine all possible neighbors up to the given order.
    The output will be: std::vector<std::pair<originalNeighbors, manyNeigbhors>>
    the many body neigbhors can be retrieved by doing:
    for (const auto nbr : manybodyNeighborIndices)
    {
        std::vector<std::pair<int,Vector3d>> neighbors = nbr.first; // this are the first orignal neighbors
        for(const auto manynbr : nbr.second)
        {
            manybody_neigbhor = neighbors;
            manybody_neigbhor.append(manynbr);
        }
    }

    this means that if originalNeigbhors.size() == 2 then for each lattice site in manyNeigbhors 
    you can combine it with originalNeigbhors to get all triplets that have these first two original neighbors (lattice indices)


    saveBothWays : bool
    if true then both i,j,k and j,i,k etc.. will be saved
    otherwise only i,j,k will be saved if i < j < k
*/

std::vector<std::pair<std::vector<std::pair<int, Vector3d>>, std::vector<std::pair<int, Vector3d>>>> ManybodyNeighborlist::build(const Neighborlist &nl, int index, int maxOrder, bool saveBothWays)
{

    std::vector<std::pair<std::vector<std::pair<int, Vector3d>>, std::vector<std::pair<int, Vector3d>>>> manybodyNeighborIndices;
    auto Ni = nl.getNeighbors(index);
    int numberOfSites = nl.size();

    int c = 2;
    Vector3d zeroVector = {0.0, 0.0, 0.0};
    std::vector<std::pair<int, Vector3d>> currentOriginalNeighbors;
    currentOriginalNeighbors.push_back(std::make_pair(index, zeroVector)); // index is always first index

    combineToHigherOrder(nl, manybodyNeighborIndices, Ni, currentOriginalNeighbors, c, saveBothWays, maxOrder);
    // }

    return manybodyNeighborIndices;
}

void ManybodyNeighborlist::combineToHigherOrder(const Neighborlist &nl,
                                                std::vector<std::pair<std::vector<std::pair<int, Vector3d>>, std::vector<std::pair<int, Vector3d>>>> &manybodyNeighborIndices,
                                                const std::vector<std::pair<int, Vector3d>> &Ni, std::vector<std::pair<int, Vector3d>> currentOriginalNeighbors, int order, bool saveBothWays, const int maxOrder)
{

    for (const auto j : Ni)
    {
        auto originalNeighborCopy = currentOriginalNeighbors;
        NeighborPairCompare comp;

        //if j is smaller than last added site then continue
        // if bothways = True then don't compare to first
        if (comp(j, originalNeighborCopy.back()) and !(saveBothWays && originalNeighborCopy.size() == 1))
        {
            continue;
        }

        originalNeighborCopy.push_back(j); // put j in originalNeigbhors

        auto N_j_offset = translateAllNi(nl.getNeighbors(j.first), j.second);
        //exclude smaller neighbors
        std::vector<std::pair<int, Vector3d>> N_j_filtered;
        for (const auto &nbrPair : N_j_offset)
        {
            if (comp(j, nbrPair)) // push back the neighbors greater than j
            {
                N_j_filtered.push_back(nbrPair);
            }
        }
        N_j_offset = N_j_filtered;

        auto intersection_N_ij = getIntersection(Ni, N_j_offset);

        if (originalNeighborCopy.size() + 1 < maxOrder)
        {
            combineToHigherOrder(nl, manybodyNeighborIndices, intersection_N_ij, originalNeighborCopy, order++, saveBothWays, maxOrder);
        }
        if (intersection_N_ij.size() > 0)
        {
            manybodyNeighborIndices.push_back(std::make_pair(originalNeighborCopy, intersection_N_ij));
        }
    }
}

std::vector<std::pair<int, Vector3d>> ManybodyNeighborlist::translateAllNi(std::vector<std::pair<int, Vector3d>> Ni, const Vector3d &unitCellOffset) const
{
    for (auto &latNbr : Ni)
    {
        latNbr.second += unitCellOffset;
    }
    return Ni;
}
