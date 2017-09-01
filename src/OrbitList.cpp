#include "OrbitList.hpp"

/**
@TODO: Think about adding a string tag here to keep track of different orbitlists
*/
OrbitList::OrbitList()
{
    //Empty constructor
}

///Construct orbitlist from mbnl and structure
OrbitList::OrbitList(const std::vector<Neighborlist> &neighborlists, const Structure &structure)
{
    _primitiveStructure = structure;
    std::unordered_map<Cluster, int> clusterIndexMap;
    ManybodyNeighborlist mbnl = ManybodyNeighborlist();

    for (int index = 0; index < structure.size(); index++)
    {
        mbnl.build(neighborlists, index, false); //bothways=false
        for (size_t i = 0; i < mbnl.getNumberOfSites(); i++)
        {
            //special case for singlet
            if (mbnl.getNumberOfSites(i) == 0)
            {
                std::vector<LatticeNeighbor> sites = mbnl.getSites(i, 0);
                Cluster cluster = Cluster(structure, sites);
                addClusterToOrbitlist(cluster, sites, clusterIndexMap);
            }

            for (size_t j = 0; j < mbnl.getNumberOfSites(i); j++)
            {
                std::vector<LatticeNeighbor> sites = mbnl.getSites(i, j);                
                Cluster cluster = Cluster(structure, sites);
                addClusterToOrbitlist(cluster, sites, clusterIndexMap);
            }
        }
    }
    bool debug = true;

    if (debug)
    {
        checkEquivalentClusters(structure);
        // std::cout << "Done checking equivalent structures" << std::endl;
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

OrbitList::OrbitList(const Structure &structure, const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, const std::vector<Neighborlist> &neighborlists)
{
    _primitiveStructure = structure;
    std::vector<std::vector<std::vector<LatticeNeighbor>>> lattice_neighbors;
    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> manybodyNeighborIndices;
    bool saveBothWays = false;
    ManybodyNeighborlist mbnl = ManybodyNeighborlist();

    //if [0,1,2] exists in taken_rows then these three rows (with columns) have been accounted for and should not be looked at
    std::unordered_set<std::vector<int>, VectorHash> taken_rows;
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

                std::vector<std::vector<LatticeNeighbor>> translatedSites = getSitesTranslatedToUnitcell(lat_nbrs);
                int missedSites = 0;

                auto sites_index_pair = getMatchesInPM(translatedSites, col1);
                auto find = taken_rows.find(sites_index_pair[0].second);
                if (find == taken_rows.end())
                {
                    //new stuff found
                    addPermutationMatrixColumns(lattice_neighbors, taken_rows, sites_index_pair[0].first, sites_index_pair[0].second, permutation_matrix, col1, true);
                }
            }

            //special singlet case
            //copy-paste from above section but with line with lat_nbrs.push_back(latnbr); is removed
            if (mbnl_pair.second.size() == 0)
            {
                std::vector<LatticeNeighbor> lat_nbrs = mbnl_pair.first;
                auto pm_rows = findRowsFromCol1(col1, lat_nbrs);
                auto find = taken_rows.find(pm_rows);
                if (find == taken_rows.end())
                {
                    //new stuff found
                    addPermutationMatrixColumns(lattice_neighbors, taken_rows, lat_nbrs, pm_rows, permutation_matrix, col1, true);
                }
            }
        }
    }

    for (int i = 0; i < lattice_neighbors.size(); i++)
    {
        std::sort(lattice_neighbors[i].begin(), lattice_neighbors[i].end());
    }

    addOrbitsFromPM(structure, lattice_neighbors);

    bool debug = true;

    if (debug)
    {
        checkEquivalentClusters(structure);
        // std::cout << "Done checking equivalent structures" << std::endl;
    }
}

/**
This will take the latticeneigbhors, and for each site outside the unitcell will translate it inside the unitcell
and translate the other sites with the same translation.

This translation will give rise to equivalent sites that sometimes are not found by using the set of crystal symmetries given
by spglib

*/
std::vector<std::vector<LatticeNeighbor>> OrbitList::getSitesTranslatedToUnitcell(const std::vector<LatticeNeighbor> &latticeNeighbors) const
{
    std::vector<std::vector<LatticeNeighbor>> translatedLatticeNeighbors;
    translatedLatticeNeighbors.push_back(latticeNeighbors);
    Vector3d zeroVector = {0.0, 0.0, 0.0};
    for (int i = 0; i < latticeNeighbors.size(); i++)
    {
        if ((latticeNeighbors[i].unitcellOffset - zeroVector).norm() > 0.1)
        {
            auto translatedSites = translateSites(latticeNeighbors, i);
            std::sort(translatedSites.begin(), translatedSites.end());

            translatedLatticeNeighbors.push_back(translatedSites);
        }
    }

    return translatedLatticeNeighbors;
}

///Take all lattice neighbors in vector latticeNeighbors and subtract the unitcelloffset of site latticeNeighbors[index]
std::vector<LatticeNeighbor> OrbitList::translateSites(const std::vector<LatticeNeighbor> &latticeNeighbors, const unsigned int index) const
{
    Vector3d offset = latticeNeighbors[index].unitcellOffset;
    auto translatedNeighbors = latticeNeighbors;
    for (auto &latNbr : translatedNeighbors)
    {
        latNbr.unitcellOffset -= offset;
    }
    return translatedNeighbors;
}

///Debug function to check that all equivalent sites in every orbit give same sorted cluster
void OrbitList::checkEquivalentClusters(const Structure &structure) const
{

    for (const auto &orbit : _orbitList)
    {
        Cluster representative_cluster = orbit.getRepresentativeCluster();
        for (const auto &sites : orbit.getEquivalentSites())
        {
            Cluster equivalentCluster = Cluster(structure, sites);
            if (representative_cluster != equivalentCluster)
            {
                std::cout << " found a \"equivalent\" cluster that were not equal representative cluster" << std::endl;
                std::cout << "representative_cluster:" << std::endl;
                representative_cluster.print();

                std::cout << "equivalentCluster:" << std::endl;
                equivalentCluster.print();

                throw std::runtime_error("found a \"equivalent\" cluster that were not equal representative cluster");
            }
        }
    }
}

/**
This adds the lattice_neighbors container found in the constructor to the orbits

each outer vector is an orbit and the inner vectors are identical sites

*/

void OrbitList::addOrbitsFromPM(const Structure &structure, const std::vector<std::vector<std::vector<LatticeNeighbor>>> &lattice_neighbors)
{

    for (const auto &equivalent_sites : lattice_neighbors)
    {
        addOrbitFromPM(structure, equivalent_sites);
    }
}

///add these equivalent sites as an orbit to orbitlist
void OrbitList::addOrbitFromPM(const Structure &structure, const std::vector<std::vector<LatticeNeighbor>> &equivalent_sites)
{

    Cluster representativeCluster = Cluster(structure, equivalent_sites[0]);
    Orbit newOrbit = Orbit(representativeCluster);
    _orbitList.push_back(newOrbit);

    for (const auto &sites : equivalent_sites)
    {
        _orbitList.back().addEquivalentSites(sites);
    }
}
/**
    From all columns in permutation matrix add all the vector<LatticeNeighbors> from pm_rows

    When taking new columns update taken_rows accordingly
 */
void OrbitList::addPermutationMatrixColumns(
    std::vector<std::vector<std::vector<LatticeNeighbor>>> &lattice_neighbors, std::unordered_set<std::vector<int>, VectorHash> &taken_rows, const std::vector<LatticeNeighbor> &lat_nbrs, const std::vector<int> &pm_rows,
    const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, const std::vector<LatticeNeighbor> &col1, bool add) const
{

    std::vector<std::vector<LatticeNeighbor>> columnLatticeNeighbors;
    columnLatticeNeighbors.reserve(permutation_matrix[0].size());
    for (size_t column = 0; column < permutation_matrix[0].size(); column++) //start at one since we got the zeroth one allready
    {
        std::vector<LatticeNeighbor> indistinctLatNbrs;

        for (const int &row : pm_rows)
        {
            indistinctLatNbrs.push_back(permutation_matrix[row][column]);
        }
        auto translatedEquivalentSites = getSitesTranslatedToUnitcell(indistinctLatNbrs);

        auto sites_index_pair = getMatchesInPM(translatedEquivalentSites, col1);

        // for (int i = 1; i < sites_index_pair.size(); i++)
        // {
        //     auto find = taken_rows.find(sites_index_pair[i].second);
        //     if( find == taken_rows.end())
        //     {

        //     }

        // }

        auto find = taken_rows.find(sites_index_pair[0].second);

        if (find == taken_rows.end())
        {
            taken_rows.insert(sites_index_pair[0].second);
            if (add)
            {
                columnLatticeNeighbors.push_back(sites_index_pair[0].first);
            }
            for (int i = 1; i < sites_index_pair.size(); i++)
            {
                find = taken_rows.find(sites_index_pair[i].second);
                if (find == taken_rows.end())
                {
                    taken_rows.insert(sites_index_pair[i].second);
                }
            }
        }
    }
    if (columnLatticeNeighbors.size() > 0)
    {
        lattice_neighbors.push_back(columnLatticeNeighbors);
    }
}

///returns the first set of translated sites that exists in col1 of permutationmatrix
std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<int>>> OrbitList::getMatchesInPM(const std::vector<std::vector<LatticeNeighbor>> &translatedSites, const std::vector<LatticeNeighbor> &col1) const
{
    std::vector<int> perm_matrix_rows;
    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<int>>> matchedSites;
    for (const auto &sites : translatedSites)
    {
        try
        {
            perm_matrix_rows = findRowsFromCol1(col1, sites);
        }
        catch (const std::runtime_error)
        {
            continue;
        }
        //no error here indicating we found matching rows in col1
        matchedSites.push_back(std::make_pair(sites, perm_matrix_rows));
    }
    if (matchedSites.size() > 0)
    {
        return matchedSites;
    }
    else
    {
        //we found no matching rows in permutation matrix, this should not happen so we throw an error

        //first print some debug info
        std::cout << "number of translated sites: " << translatedSites.size() << std::endl;
        std::cout << "sites: " << std::endl;
        for (auto latnbrs : translatedSites)
        {
            for (auto latnbr : latnbrs)
            {
                latnbr.print();
            }
            std::cout << " ========= " << std::endl;
        }
        std::cout << "col1:" << std::endl;
        for (auto row : col1)
        {
            row.print();
        }
        throw std::runtime_error("Did not find any of the translated sites in col1 of permutation matrix in function getFirstMatchInPM in orbitlist");
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
                //latNbrp.print();
            }
            //  latNbr.print();
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
