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







OrbitList::OrbitList(const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, const std::vector<Neighborlist> &neighborlists)
{

    std::vector<std::vector<std::vector<LatticeNeighbor>>> lattice_neighbors;
    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> manybodyNeighborIndices;
    bool saveBothWays = false;
    ManybodyNeighborlist mbnl = ManybodyNeighborlist();

    //if [0,1,2] exists in taken_rows then these three rows (with columns) have been accounted for and should not be looked at
    std::vector<std::vector<int>> taken_rows;
    std::vector<LatticeNeighbor> col1 = getColumn1FromPM(permutation_matrix, false);

    std::set<LatticeNeighbor> col1_uniques(col1.begin(), col1.end());
    if (col1.size() != col1_uniques.size())
    {
        throw std::runtime_error("Found duplicates in column1 of permutation matrix");
    }

    std::cout << "Size of col1 " << col1.size() << std::endl;
    for (size_t index = 0; index < neighborlists[0].size(); index++)
    {

        std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> mbnl_latnbrs = mbnl.build(neighborlists, index, saveBothWays);
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

 
}

/**
    From all columns in permutation matrix add all the vector<LatticeNeighbors> from pm_rows

    When taking new columns update taken_rows accordingly
 */
void OrbitList::addPermutationMatrixColumns(
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

bool OrbitList::validatedCluster(const std::vector<LatticeNeighbor> &latticeNeighbors) const
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
std::vector<int> OrbitList::findRowsFromCol1(const std::vector<LatticeNeighbor> &col1, const std::vector<LatticeNeighbor> &latNbrs, bool sortIt) const
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
std::vector<LatticeNeighbor> OrbitList::getColumn1FromPM(const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, bool sortIt) const
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
