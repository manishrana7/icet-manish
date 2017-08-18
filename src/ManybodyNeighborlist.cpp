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

std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> ManybodyNeighborlist::build(const std::vector<Neighborlist> &neighborlists, int index, bool saveBothWays)
{
    
    if(neighborlists.empty())
    {
        throw std::runtime_error("Error: neigbhorlist vector is empty in ManybodyNeighborlist::build");
    }
    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> manybodyNeighborIndices;

    addSinglet(index, manybodyNeighborIndices);
    addPairs(index, neighborlists[0], manybodyNeighborIndices, saveBothWays);

    for (size_t c = 2; c < neighborlists.size() + 2; c++)
    {
        auto Ni = neighborlists[c - 2].getNeighbors(index);
        int numberOfSites = neighborlists[c - 2].size();
        Vector3d zeroVector = {0.0, 0.0, 0.0};
        std::vector<LatticeNeighbor> currentOriginalNeighbors;
        currentOriginalNeighbors.push_back(LatticeNeighbor(index, zeroVector)); // index is always first index

        combineToHigherOrder(neighborlists[c - 2], manybodyNeighborIndices, Ni, currentOriginalNeighbors, saveBothWays, c);
    }
    return manybodyNeighborIndices;
}

///Adds singlet from the index to manybodyNeighborIndices
void ManybodyNeighborlist::addSinglet(const int index, std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &manybodyNeighborIndices) const
{
    Vector3d zeroVector = {0.0, 0.0, 0.0};
    LatticeNeighbor latticeNeighborSinglet = LatticeNeighbor(index, zeroVector);
    std::vector<LatticeNeighbor> singletLatticeNeighbors;
    singletLatticeNeighbors.push_back(latticeNeighborSinglet);

    std::vector<LatticeNeighbor> latNbrsEmpty;
    manybodyNeighborIndices.push_back(std::make_pair(singletLatticeNeighbors, latNbrsEmpty));
}

///Add all pairs originating from index using neighborlist
void ManybodyNeighborlist::addPairs(const int index, const Neighborlist &neighborList,
                                    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &manybodyNeighborIndices, bool saveBothWays) const

{
    Vector3d zeroVector = {0.0, 0.0, 0.0};
    LatticeNeighbor latticeNeighborIndex = LatticeNeighbor(index, zeroVector);
    
    std::vector<LatticeNeighbor> firstSite = {latticeNeighborIndex};
    std::vector<LatticeNeighbor> Ni = neighborList.getNeighbors(index);
    //exclude smaller neighbors
    if (!saveBothWays)
    {
        Ni = getFilteredNj(Ni, latticeNeighborIndex);
    }

    if(Ni.size()==0)
    {
        return;
    }
    manybodyNeighborIndices.push_back(std::make_pair(firstSite, Ni));
}

/**
This will use the permutationmatrix together with the neighborlist to construct the distinct and indistinct sets of points

The output will be std::vector<std::vector<std::vector<LatticeNeighbor>>> 
the next outer vector contains the set of indistinct set of lattice neighbors.


Algorithm
=========

The algorithm will work by taking the first column of permutation matrix: col1 = permutation_matrix[:,0]

and it will take Ni (lattice neighbors of i) and find the intersection of Ni and col1, intersection(Ni, col1) = Ni_pm
all j in Ni_c1 are then within the cutoff of site i, then depending on the order all the pairs/triplets will be constructed from the 
lattice neighbors in Ni_pm.

This will be repeated for all the sites in the neighborlist. In the end all the pair/triplet terms will have been generated from col1.

Then you will take the first vector<LatticeNeighbors> and find the rows of these LatNbrs in col1 of the permutation matrix.
then you traverse through all the columns in the permutation matrix. These vector<latticeNeighbors> are then indistinct from the original. 
Note that here a validity check is made to ensure that atleast one LatticeNeigbhor originate in the original lattice (unitcell offset = [0,0,0])
otherwise we are including a  "ghost cluster".

The new vector<latticeNeighbors> found when traversing the columns are likely to have been found from the combinations in Ni_pm and these must then 
be removed/overlooked when moving to the next vector<LatticeNeighbor>.

*/

std::vector<std::vector<std::vector<LatticeNeighbor>>> ManybodyNeighborlist::buildFromPermutationMatrix(const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, const std::vector<Neighborlist> &neighborlists)
{

    std::vector<std::vector<std::vector<LatticeNeighbor>>> lattice_neighbors;
    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> manybodyNeighborIndices;
    bool saveBothWays = false;

    //if [0,1,2] exists in taken_rows then these three rows (with columns) have been accounted for and should not be looked at
    std::vector<std::vector<int>> taken_rows;
    std::vector<LatticeNeighbor> col1 = getColumn1FromPM(permutation_matrix, false);

    std::set<LatticeNeighbor> col1_uniques(col1.begin(), col1.end());
    if (col1.size() != col1_uniques.size())
    {
        throw std::runtime_error("Found duplicates in column1 of permutation matrix");
    }
    // for (auto latnbr : col1_uniques)
    // {
    //     latnbr.print();
    // }
    std::cout << "Size of col1 " << col1.size() << std::endl;
    for (size_t index = 0; index < neighborlists[0].size(); index++)
    {

        std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> mbnl_latnbrs = build(neighborlists, index, saveBothWays);
        for (const auto &mbnl_pair : mbnl_latnbrs)
        {

            for (const auto &latnbr : mbnl_pair.second)
            {
                std::vector<LatticeNeighbor> lat_nbrs = mbnl_pair.first;
                lat_nbrs.push_back(latnbr);
                // std::cout<<"find in build function:"<<std::endl;
                auto pm_rows = findRowsFromCol1(col1, lat_nbrs);
                auto find = std::find(taken_rows.begin(), taken_rows.end(), pm_rows);
                if (find == taken_rows.end())
                {
                    //new stuff found
                    addPermutationMatrixColumns(lattice_neighbors, taken_rows, lat_nbrs, pm_rows, permutation_matrix, col1);
                }
            }
        }
    }

    for (int i = 0; i < lattice_neighbors.size(); i++)
    {
        std::sort(lattice_neighbors[i].begin(), lattice_neighbors[i].end());
    }

    for (int i = 0; i < lattice_neighbors.size(); i++)
    {
        for (int j = i + 1; j < lattice_neighbors.size(); j++)
        {
            if (lattice_neighbors[i] == lattice_neighbors[j])
            {
                std::cout << "Found duplicate" << std::endl;
            }
        }
    }
    std::sort(taken_rows.begin(), taken_rows.end());
    std::cout << "Taken rows: " << std::endl;
    for (auto taken_row : taken_rows)
    {
        for (auto el : taken_row)
        {
            std::cout << el << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "=========" << std::endl;

    return lattice_neighbors;
}

/**
    From all columns in permutation matrix add all the vector<LatticeNeighbors> from pm_rows

    When taking new columns update taken_rows accordingly
 */
void ManybodyNeighborlist::addPermutationMatrixColumns(
    std::vector<std::vector<std::vector<LatticeNeighbor>>> &lattice_neighbors, std::vector<std::vector<int>> &taken_rows, const std::vector<LatticeNeighbor> &lat_nbrs, const std::vector<int> &pm_rows,
    const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, const std::vector<LatticeNeighbor> &col1) const
{

    std::vector<std::vector<LatticeNeighbor>> columnLatticeNeighbors;
    columnLatticeNeighbors.reserve(permutation_matrix[0].size());
    columnLatticeNeighbors.push_back(lat_nbrs);

    taken_rows.push_back(pm_rows);
    for (size_t column = 1; column < permutation_matrix[0].size(); column++) //start at one since we got the zeroth one allready
    {
        std::vector<LatticeNeighbor> indistinctLatNbrs;

        for (const int &row : pm_rows)
        {
            indistinctLatNbrs.push_back(permutation_matrix[row][column]);
        }
        auto perm_matrix_rows = findRowsFromCol1(col1, indistinctLatNbrs);
        // std::cout<<"find in permutation matrix function:"<<std::endl;
        auto find = std::find(taken_rows.begin(), taken_rows.end(), perm_matrix_rows);
        if (find == taken_rows.end())
        {
            // std::cout<<"Found new taken rows: "<<std::endl;
            // for(auto el : perm_matrix_rows){std::cout<<el<<" ";}
            // std::cout<<std::endl;
            // std::cout<<"Taken rows: "<<std::endl;
            // for(auto taken_row : taken_rows)
            // {
            //     for(auto el : taken_row){std::cout<<el<<" ";}
            //     std::cout<<std::endl;
            // }
            // std::cout<<"========="<<std::endl;
            taken_rows.push_back(perm_matrix_rows);
        }
        columnLatticeNeighbors.push_back(indistinctLatNbrs);
    }
    if (columnLatticeNeighbors.size() > 0)
    {
        lattice_neighbors.push_back(columnLatticeNeighbors);
    }
}

/**

Checks that atleast one lattice neigbhor originate in the original cell (has one cell offset = [0,0,0])
*/

bool ManybodyNeighborlist::validatedCluster(const std::vector<LatticeNeighbor> &latticeNeighbors) const
{
    Vector3d zeroVector = {0., 0., 0.};
    for (const auto &latNbr : latticeNeighbors)
    {
        if (latNbr.unitcellOffset == zeroVector)
        {
            return true;
        }
    }
    return false;
}

/**
 Searches for latticeNeighbors in col1 of permutation matrix and find the corresponding rows 
*/
std::vector<int> ManybodyNeighborlist::findRowsFromCol1(const std::vector<LatticeNeighbor> &col1, const std::vector<LatticeNeighbor> &latNbrs, bool sortIt) const
{
    std::vector<int> rows;
    for (const auto &latNbr : latNbrs)
    {
        const auto find = std::find(col1.begin(), col1.end(), latNbr);
        if (find == col1.end())
        {
            for (const auto &latNbrp : latNbrs)
            {
                latNbrp.print();
            }
            latNbr.print();
            throw std::runtime_error("Did not find lattice neigbhor in col1 of permutation matrix in function findRowsFromCol1 in mbnl");
        }
        else
        {
            int row_in_col1 = std::distance(col1.begin(), find);
            rows.push_back(row_in_col1);
        }
    }
    // if (sortIt)
    // {
    std::sort(rows.begin(), rows.end());
    // }
    return rows;
}

/**
    Returns the first column of the permutation matrix

    named arguments:
        sortIt : if true it will sort col1 (default true)

*/
std::vector<LatticeNeighbor> ManybodyNeighborlist::getColumn1FromPM(const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, bool sortIt) const
{
    std::vector<LatticeNeighbor> col1;
    col1.reserve(permutation_matrix[0].size());
    for (const auto &row : permutation_matrix)
    {
        col1.push_back(row[0]);
    }
    if (sortIt)
    {
        std::sort(col1.begin(), col1.end());
    }
    return col1;
}

/*
    for each j in Ni construct the intersect of N_j and N_i = N_ij.
    all k in N_ij are then neighbors with i,j
    what is saved is then i,j and N_ij up to the desired order "maxorder"
*/
void ManybodyNeighborlist::combineToHigherOrder(const Neighborlist &nl,
                                                std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &manybodyNeighborIndices,
                                                const std::vector<LatticeNeighbor> &Ni, std::vector<LatticeNeighbor> &currentOriginalNeighbors, bool saveBothWays, const int maxOrder)
{

    for (const auto &j : Ni)
    {
        //if j is smaller than last added site then continue
        // if bothways = True then don't compare to first
        bool cont = false;

        if (saveBothWays)
        {
            if (currentOriginalNeighbors.size() > 1)
            {
                if (j < currentOriginalNeighbors.back())
                {
                    cont = true;
                }
            }
        }
        else
        {
            if (j < currentOriginalNeighbors.back())
            {
                cont = true;
            }
        }
        if (cont)
        {
            continue;
        }
        // if ((!saveBothWays && currentOriginalNeighbors.size() == 1) && j < currentOriginalNeighbors.back())
        // {
        //     continue;
        // }

        auto originalNeighborCopy = currentOriginalNeighbors;
        originalNeighborCopy.push_back(j); // put j in originalNeigbhors

        auto Nj = nl.getNeighbors(j.index);

        //translate the neighbors
        translateAllNi(Nj, j.unitcellOffset);

        //exclude smaller neighbors
        if (!saveBothWays)
        {
            Nj = getFilteredNj(Nj, j);
        }

        //construct the intersection
        const auto intersection_N_ij = getIntersection(Ni, Nj);

        if (originalNeighborCopy.size() + 1 < maxOrder)
        {
            combineToHigherOrder(nl, manybodyNeighborIndices, intersection_N_ij, originalNeighborCopy, saveBothWays, maxOrder);
        }

        if (intersection_N_ij.size() > 0 && originalNeighborCopy.size() == (maxOrder - 1))
        {
            manybodyNeighborIndices.push_back(std::make_pair(originalNeighborCopy, intersection_N_ij));
        }
    }
}

/*
Since N_j is always sorted then simply search for first k in N_j that have k>= j
and then filtered are from indexof(k) to end()

*/
std::vector<LatticeNeighbor> ManybodyNeighborlist::getFilteredNj(const std::vector<LatticeNeighbor> &N_j, const LatticeNeighbor &j) const
{
    auto first = std::upper_bound(N_j.begin(), N_j.end(), j);

    std::vector<LatticeNeighbor> ret(first, N_j.end());
    return ret;
}

/**
    Offsets all indice, offsets pairs in Ni with the input offset, e.g:
    For all j in Ni:
     offset j.offset with "unitCellOffset"

*/
void ManybodyNeighborlist::translateAllNi(std::vector<LatticeNeighbor> &Ni, const Vector3d &unitCellOffset) const
{
    for (auto &latNbr : Ni)
    {
        latNbr.unitcellOffset += unitCellOffset;
    }
}
