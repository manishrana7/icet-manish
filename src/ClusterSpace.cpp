#include "ClusterSpace.hpp"

//namespace icet {

using namespace std;

/**
@details Initializes a ClusterSpace object.
@param numberOfAllowedComponents number of allowed components for each site of the primitive structure
@param chemicalSymbols chemical symbol for each site
@param orbitList list of orbits for the primitive structure
*/
ClusterSpace::ClusterSpace(vector<int> numberOfAllowedComponents,
                           vector<string> chemicalSymbols,
                           const OrbitList orbitList)
{
    _numberOfAllowedComponents = numberOfAllowedComponents;
    _orbitList = orbitList;
    _primitiveStructure = orbitList.getPrimitiveStructure();
    _primitiveStructure.setNumberOfAllowedComponents(_numberOfAllowedComponents);
    /// @todo Why is the task of setting up the element map outsourced to a private function that is not used anywhere else?
    setupElementMap(chemicalSymbols);
    /// @todo Why is the collectClusterSpaceInfo function not executed immediately, which would render _isClusterSpaceInitialized unnecessary?
    _isClusterSpaceInitialized = false;
};


/**
@brief Sets up a map between chemical elements and the internal species enumeration scheme.
@param elements list of chemical symbols
*/
void ClusterSpace::setupElementMap(const vector<string> elements)
{
    // convert chemical symbols to atomic numbers
    vector<int> intElements;
    for (const auto el : elements)
        intElements.push_back(PeriodicTable::strInt[el]);

    // sort list of elements
    sort(intElements.begin(), intElements.end());
    _elements = intElements;

    // map atomic numbers to internal species enumeration scheme
    for (size_t i = 0; i < elements.size(); i++)
        _elementMap[intElements[i]] = i;
}


/**
@details Calculates and then returns the cluster vector for the input
structure in the cluster space. The first element in the cluster vector will
always be one (1) corresponding to the zerolet. The remaining elements of the
cluster vector represent averages over orbits (symmetry equivalent clusters)
of increasing order and size.

@param structure input configuration

@todo review the necessity for having the orderIntact argument to e.g., countOrbitList.
**/
vector<double> ClusterSpace::getClusterVector(const Structure &structure) const
{
    // count the clusters in the orbit with the same orientation (order) as the prototype cluster
    // @todo Clarify description.
    bool orderIntact = true;

    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_orbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();

    // Create local orbit list and count associated clusters.
    for (int i = 0; i < uniqueOffsets; i++) {
        const auto localOrbitList = localOrbitListGenerator.getLocalOrbitList(i);
        clusterCounts.countOrbitList(structure, localOrbitList, orderIntact);
    }

    // Check that the number of unique offsets equals the number of unit cells in the supercell.
    int numberOfUnitcellRepetitions = structure.size() / _primitiveStructure.size();
    if (uniqueOffsets != numberOfUnitcellRepetitions) {
        string errorMessage = "The number of unique offsets does not match the number of primitive units in the input structure: ";
        errorMessage += to_string(uniqueOffsets) + " != " + to_string(numberOfUnitcellRepetitions);
        throw runtime_error(errorMessage);
    }

    /// @todo Describe what is happening here.
    const auto clusterMap = clusterCounts.getClusterCounts();
    vector<double> clusterVector;

    // insert the zerolet.
    clusterVector.push_back(1);

// CONTINUE HERE

    // Loop over orbits.
    /// @todo Turn this into a proper loop over orbits.
    for (size_t i = 0; i < _orbitList.size(); i++) {

        auto repCluster = _orbitList.getOrbit(i).getRepresentativeCluster();
        vector<int> numberOfAllowedOccupations;
        try {
            numberOfAllowedOccupations = getNumberOfAllowedComponentsForEachSite(_primitiveStructure, _orbitList.getOrbit(i).getRepresentativeSites());
        }
        catch (const exception& e) {
            throw runtime_error("Failed retrieving the number of allowed occupations in ClusterSpace::getClusterVector");
        }

        // Jump to the next orbit if none of the sites in the representative cluster are active (i.e. there allowed occupation is less than 2).
        if (any_of(numberOfAllowedOccupations.begin(), numberOfAllowedOccupations.end(), [](int numberOfAllowedOccupation){ return numberOfAllowedOccupation < 2; }))
            continue;

        auto mcVectors = _orbitList.getOrbit(i).getMCVectors(numberOfAllowedOccupations);
        auto allowedPermutationsSet = _orbitList.getOrbit(i).getAllowedSitesPermutations();
        auto elementPermutations = getMultiComponentVectorPermutations(mcVectors, i);
        repCluster.setClusterTag(i);
        int currentMCVectorIndex = 0;
        for (const auto &mcVector : mcVectors)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            for (const auto &elementsCountPair : clusterMap.at(repCluster))
            {

                // TODO check if numberOfAllowedOccupations should be permuted as well.
                for (const auto &perm : elementPermutations[currentMCVectorIndex])
                {
                    auto permutedMCVector = icet::getPermutedVector(mcVector, perm);
                    auto permutednumberOfAllowedOccupations = icet::getPermutedVector(numberOfAllowedOccupations, perm);
                    clusterVectorElement += getClusterProduct(permutedMCVector, permutednumberOfAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            clusterVectorElement /= ((double)multiplicity);
            clusterVector.push_back(clusterVectorElement);

            currentMCVectorIndex++;
        }
    }
    return clusterVector;
}


/// Returns the native clusters count in this structure, i.e. only clusters inside the unit cell
ClusterCounts ClusterSpace::getNativeClusters(const Structure &structure) const
{
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_orbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();
    int tags = 0;
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto localOrbitList = localOrbitListGenerator.getLocalOrbitList(i);
        for (int j = 0; j < localOrbitList.size(); j++)
        {
            Cluster repr_cluster = localOrbitList.getOrbit(j).getRepresentativeCluster();

            for (const auto sites : localOrbitList.getOrbit(j).getPermutedEquivalentSites())
            {
                bool cont = true;
                for (auto site : sites)
                {
                    if (site.unitcellOffset().norm() > 0.1)
                    {
                        cont = false;
                    }
                }
                if (!cont)
                {
                    continue;
                }
                repr_cluster.setClusterTag(j);
                if (repr_cluster.order() != 1)
                {
                    vector<int> elements(sites.size());
                    for (size_t i = 0; i < sites.size(); i++)
                    {
                        elements[i] = structure.getAtomicNumber(sites[i].index());
                    }
                    clusterCounts.countCluster(repr_cluster, elements, orderIntact);
                }
                else
                {
                    vector<int> elements(sites.size());
                    for (size_t i = 0; i < sites.size(); i++)
                    {
                        elements[i] = structure.getAtomicNumber(sites[i].index());
                    }
                    clusterCounts.countCluster(repr_cluster, elements, orderIntact);
                }
            }
        }
    }
    return clusterCounts;
}

/**
  @details This method return the mc vector permutations for each mc vector.
  Example1: Given mc vectors [0, 0], [0,1] and [1,1]
  the returned permutations should be [[1,0]], [[0,1],[1,0]], [1,1].
  i.e. the [0,1] mc vector should count elements with permutations [1,0] and [1,0]

  Given mc vectors [0, 0], [0,1], [1,0] and [1,1] the returned permutations
  will only be the self permutations since the mc vectors [0,1] and [1,0] will handle
  the AB vs BA choice.

  @param mcVectors the mc vectors for this orbit
  @param orbitIndex : The orbit index to take the allowed permutations from.

*/

vector<vector<vector<int>>> ClusterSpace::getMultiComponentVectorPermutations(const vector<vector<int>> &mcVectors, const int orbitIndex) const
{
    const auto allowedPermutations = _orbitList.getOrbit(orbitIndex).getAllowedSitesPermutations();

    vector<vector<vector<int>>> elementPermutations;
    vector<int> selfPermutation;
    for (int i = 0; i < mcVectors[0].size(); i++)
    {
        selfPermutation.push_back(i);
    }

    for (const auto &mc : mcVectors)
    {
        vector<vector<int>> mcPermutations;
        mcPermutations.push_back(selfPermutation);
        vector<vector<int>> takenPermutations;
        takenPermutations.push_back(selfPermutation);
        for (const vector<int> perm : allowedPermutations)
        {
            auto permutedMcVector = icet::getPermutedVector(mc, perm);
            auto findPerm = find(mcVectors.begin(), mcVectors.end(), permutedMcVector);
            auto findIfTaken = find(takenPermutations.begin(), takenPermutations.end(), permutedMcVector);
            if (findPerm == mcVectors.end() && findIfTaken == takenPermutations.end() && mc != permutedMcVector)
            {
                mcPermutations.push_back(perm);
                takenPermutations.push_back(permutedMcVector);
            }
        }
        sort(mcPermutations.begin(), mcPermutations.end());
        elementPermutations.push_back(mcPermutations);
    }
    return elementPermutations;
}

/**
This is the default cluster function
*/
double ClusterSpace::defaultClusterFunction(const int numberOfAllowedComponents, const int clusterFunction, const int element) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)numberOfAllowedComponents));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)element / ((double)numberOfAllowedComponents));
    }
}

/// Returns the full cluster product of the entire cluster (elements vector) assuming all sites have same numberOfAllowedComponents.
double ClusterSpace::getClusterProduct(const vector<int> &mcVector, const vector<int> &numberOfAllowedComponents, const vector<int> &elements) const
{
    double clusterProduct = 1;
    for (int i = 0; i < elements.size(); i++)
    {
        clusterProduct *= defaultClusterFunction(numberOfAllowedComponents[i], mcVector[i], _elementMap.at(elements[i]));
    }
    return clusterProduct;
}

/// Returns the allowed occupations on the sites
/// @todo Why is function not a member of Structure?
vector<int> ClusterSpace::getNumberOfAllowedComponentsForEachSite(const Structure &structure, const vector<LatticeSite> &latticeNeighbors) const
{
    vector<int> numberOfAllowedComponents;
    numberOfAllowedComponents.reserve(latticeNeighbors.size());
    for (const auto &latnbr : latticeNeighbors)
    {
        numberOfAllowedComponents.push_back(structure.getNumberOfAllowedComponentsForEachSite(latnbr.index()));
    }
    return numberOfAllowedComponents;
}

/// Collect information about the cluster space.
void ClusterSpace::collectClusterSpaceInfo()
{
    _clusterSpaceInfo.clear();
    vector<int> emptyVec = {0};
    _clusterSpaceInfo.push_back(make_pair(-1, emptyVec));
    for (int i = 0; i < _orbitList.size(); i++)
    {
        auto numberOfAllowedOccupations = getNumberOfAllowedComponentsForEachSite(_primitiveStructure, _orbitList.getOrbit(i).getRepresentativeSites());
        auto mcVectors = _orbitList.getOrbit(i).getMCVectors(numberOfAllowedOccupations);
        for (const auto &mcVector : mcVectors)
        {
            _clusterSpaceInfo.push_back(make_pair(i, mcVector));
        }
    }
    _isClusterSpaceInitialized = true;
}

/// Returns information about the cluster space.
pair<int, vector<int>> ClusterSpace::getClusterSpaceInfo(const unsigned int index)
{
    if (!_isClusterSpaceInitialized)
    {
        collectClusterSpaceInfo();
    }

    if (index >= _clusterSpaceInfo.size())
    {
        string errMSG = "Out of range in ClusterSpace::getClusterSpaceInfo " + to_string(index) + " >= " + to_string(_clusterSpaceInfo.size());
        throw out_of_range(errMSG);
    }

    return _clusterSpaceInfo[index];
}

//}