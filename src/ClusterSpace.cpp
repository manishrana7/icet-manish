#include "ClusterSpace.hpp"

//namespace icet {

/**
@details This constructor initializes a ClusterSpace object.
@param numberOfAllowedSpecies number of allowed components for each site of the primitive structure
@param chemicalSymbols chemical symbol for each site
@param orbitList list of orbits for the primitive structure
*/
ClusterSpace::ClusterSpace(std::vector<int> numberOfAllowedSpecies,
                           std::vector<std::string> chemicalSymbols,
                           const OrbitList orbitList)
{
    _numberOfAllowedSpeciesPerSite = numberOfAllowedSpecies;
    _orbitList = orbitList;
    _primitiveStructure = orbitList.getPrimitiveStructure();
    _primitiveStructure.setNumberOfAllowedSpecies(_numberOfAllowedSpeciesPerSite);

    // Set up a map between chemical species and the internal species enumeration scheme.
    for (const auto el : chemicalSymbols)
    {
        _species.push_back(PeriodicTable::strInt[el]);
    }
    sort(_species.begin(), _species.end());
    for (size_t i = 0; i < _species.size(); i++)
    {
        _speciesMap[_species[i]] = i;
    }


    //Remove orbits that are inactive
    for(int i=_orbitList.size()-1; i>=0; i--)
    {
        auto numberOfAllowedSpecies = getNumberOfAllowedSpeciesBySite(_primitiveStructure, _orbitList.getOrbit(i).getRepresentativeSites());

        if (std::any_of(numberOfAllowedSpecies.begin(), numberOfAllowedSpecies.end(), [](int n) { return n < 2; }))
        {
            _orbitList.removeOrbit(i);
        }        
    }


    /// @todo Why is the collectClusterSpaceInfo function not executed
    /// immediately, which would render _isClusterSpaceInitialized unnecessary?
    _isClusterSpaceInitialized = false;
    collectClusterSpaceInfo();


    
}

/**
@details This function calculates and then returns the cluster vector for the
input structure in the cluster space.

The first element in the cluster vector will always be one (1) corresponding to
the zerolet. The remaining elements of the cluster vector represent averages
over orbits (symmetry equivalent clusters) of increasing order and size.

@param structure input configuration

@todo review the necessity for having the orderIntact argument to e.g.,
countOrbitList.
**/
std::vector<double> ClusterSpace::getClusterVector(const Structure &structure) const
{
    // count the clusters in the orbit with the same orientation (order) as the prototype cluster
    // @todo Clarify description.
    bool orderIntact = true;
    bool permuteSites = true;

    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(_orbitList, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getNumberOfUniqueOffsets();
    ClusterCounts clusterCounts = ClusterCounts();

    // Create local orbit list and count associated clusters.
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto localOrbitList = localOrbitListGenerator.getLocalOrbitList(i);
        clusterCounts.countOrbitList(structure, localOrbitList, orderIntact, permuteSites);
    }

    // Check that the number of unique offsets equals the number of unit cells in the supercell.
    int numberOfUnitcellRepetitions = structure.size() / _primitiveStructure.size();
    if (uniqueOffsets != numberOfUnitcellRepetitions)
    {
        std::string msg = "The number of unique offsets does not match the number of primitive units in the input structure: ";
        msg += std::to_string(uniqueOffsets) + " != " + std::to_string(numberOfUnitcellRepetitions);
        throw std::runtime_error(msg);
    }

    /// @todo Describe what is happening here.
    const auto clusterMap = clusterCounts.getClusterCounts();

    // Initialize cluster vector and insert zerolet.
    std::vector<double> clusterVector;
    clusterVector.push_back(1);

    // Loop over orbits.
    /// @todo Turn this into a proper loop over orbits. This probably requires upgrading OrbitList.
    for (size_t i = 0; i < _orbitList.size(); i++)
    {

        auto representativeCluster = _orbitList._orbitList[i].getRepresentativeCluster();
        // @todo This is necessary. Side effects need to carefully evaluated. Ideally move this task elsewhere as it is repeated for every structure, for which the cluster vector is computed.
        representativeCluster.setTag(i);

        std::vector<int> numberOfAllowedSpeciesBySite;
        try
        {
            numberOfAllowedSpeciesBySite = getNumberOfAllowedSpeciesBySite(_primitiveStructure, _orbitList._orbitList[i].getRepresentativeSites());
        }
        catch (const std::exception &e)
        {
            std::string msg = "Failed retrieving the number of allowed species in ClusterSpace::getClusterVector\n";
            msg += e.what();
            throw std::runtime_error(msg);
        }

        // Jump to the next orbit if any of the sites in the representative cluster are in-active (i.e. the number of allowed species on this site is less than 2).
        if (any_of(numberOfAllowedSpeciesBySite.begin(), numberOfAllowedSpeciesBySite.end(), [](int n) { return n < 2; }))
        {
            continue;
        }

        // First we obtain the multi-component vectors for this orbit, i.e. a vector
        // of vectors of int (where the int represents a cluster function index).
        // Example 1: For an AB alloy we obtain [0, 0] and [0, 0, 0] for pair and triplet terms, respectively.
        // Example 2: For an ABC alloy we obtain [0, 0], [0, 1], [1, 1] for pairs and similarly for triplets.
        // Depending on the symmetry of the cluster one might also obtain [1, 0] (e.g., in a clathrate or for some clusters on a HCP lattice).

        // auto multiComponentVectors = _orbitList.getOrbit(i).getMultiComponentVectors(numberOfAllowedSpeciesBySite);
        const auto &multiComponentVectors = _multiComponentVectors[i];

        // @todo Make getMultiComponentVectorPermutations take an Orbit rather than an index. Then swap the loop over int for a loop over Orbit above.
        // auto speciesPermutations = getMultiComponentVectorPermutations(multiComponentVectors, i);
        const auto &speciesPermutations = _elementPermutations[i];
        int currentMultiComponentVectorIndex = 0;
        for (const auto &multiComponentVector : multiComponentVectors)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            for (const auto &speciesCountPair : clusterMap.at(representativeCluster))
            {

                /// @todo Check if numberOfAllowedSpecies should be permuted as well. Is this todo still relevant?
                for (const auto &perm : speciesPermutations[currentMultiComponentVectorIndex])
                {
                    auto permutedMultiComponentVector = icet::getPermutedVector(multiComponentVector, perm);
                    auto permutedNumberOfAllowedSpeciesBySite = icet::getPermutedVector(numberOfAllowedSpeciesBySite, perm);
                    clusterVectorElement += evaluateClusterProduct(permutedMultiComponentVector, permutedNumberOfAllowedSpeciesBySite, speciesCountPair.first) * speciesCountPair.second;
                    multiplicity += speciesCountPair.second;
                }
            }
            clusterVectorElement /= ((double)multiplicity);
            clusterVector.push_back(clusterVectorElement);

            currentMultiComponentVectorIndex++;
        }
    }
    return clusterVector;
}

/**
  @details This method return the multi-component vector permutations for each
  multi-component vector.

  Example 1: Given multi-component vectors [0, 0], [0, 1] and [1, 1]
  the returned permutations should be [[1, 0]], [[0, 1],[1, 0]], [1, 1].
  i.e. the [0, 1] multi-component vector should count elements with
  permutations [1, 0] and [1, 0].

  Example 2: Given multi-component vectors [0, 0], [0, 1], [1, 0] and [1, 1]
  the returned permutations will only be the self permutations since the
  multi-component vectors [0, 1] and [1, 0] will handle the AB vs BA choice.

  @param multiComponentVectors multi-component vectors for this orbit
  @param orbitIndex index from which to take the allowed permutations

  @returns a vector of a vector of a vector of ints; here the innermost index

  @todo This function should take an Orbit rather than an orbit index.
*/

std::vector<std::vector<std::vector<int>>> ClusterSpace::getMultiComponentVectorPermutations(const std::vector<std::vector<int>> &multiComponentVectors, const int orbitIndex) const
{
    const auto allowedPermutations = _orbitList._orbitList[orbitIndex].getAllowedSitesPermutations();

    std::vector<std::vector<std::vector<int>>> elementPermutations;
    std::vector<int> selfPermutation;
    for (int i = 0; i < multiComponentVectors[0].size(); i++)
    {
        selfPermutation.push_back(i);
    }

    for (const auto &mc : multiComponentVectors)
    {
        std::vector<std::vector<int>> mcPermutations;
        mcPermutations.push_back(selfPermutation);
        std::vector<std::vector<int>> takenPermutations;
        takenPermutations.push_back(selfPermutation);
        for (const std::vector<int> perm : allowedPermutations)
        {
            auto permutedMultiComponentVector = icet::getPermutedVector(mc, perm);
            auto findPerm = find(multiComponentVectors.begin(), multiComponentVectors.end(), permutedMultiComponentVector);
            auto findIfTaken = find(takenPermutations.begin(), takenPermutations.end(), permutedMultiComponentVector);
            if (findPerm == multiComponentVectors.end() && findIfTaken == takenPermutations.end() && mc != permutedMultiComponentVector)
            {
                mcPermutations.push_back(perm);
                takenPermutations.push_back(permutedMultiComponentVector);
            }
        }
        sort(mcPermutations.begin(), mcPermutations.end());
        elementPermutations.push_back(mcPermutations);
    }
    return elementPermutations;
}

/**
@details Evaluates the cluster function using the specified parameters.

The cluster functions (also "orthogonal point functions") are defined as

.. math::

   \Theta_{n}(\sigma_p) = \begin{cases}
      1                                    &\quad \text{if}~n=0 \\
      -\cos\left(\pi(n+1)\sigma_p/M\right) &\quad \text{if n is odd} \\
      -\sin\left(\pi n   \sigma_p/M\right) &\quad \text{if n is even}
    \end{cases}

@param numberOfAllowedSpecies number of allowed species on the site in question
@param clusterFunction index of cluster function
@param species index of species

@returns the value of the cluster function
*/
double ClusterSpace::evaluateClusterFunction(const int numberOfAllowedSpecies, const int clusterFunction, const int species) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)species / ((double)numberOfAllowedSpecies));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) * (double)species / ((double)numberOfAllowedSpecies));
    }
}

/**
@details Evaluates the full cluster product of the entire cluster.

@param multiComponentVector multi-component vector, each element of the vector gives the index of a cluster function
@param numberOfAllowedSpecies number of species allowed on the sites in this cluster (all sites involved are assumed to have the same number of allowed species)
@param species species that occupy (decorate) the cluster identified by atomic number

@returns the cluster product
**/
double ClusterSpace::evaluateClusterProduct(const std::vector<int> &multiComponentVector, const std::vector<int> &numberOfAllowedSpecies, const std::vector<int> &species) const
{
    double clusterProduct = 1;
    for (int i = 0; i < species.size(); i++)
    {
        clusterProduct *= evaluateClusterFunction(numberOfAllowedSpecies[i], multiComponentVector[i], _speciesMap.at(species[i]));
    }
    return clusterProduct;
}

/**
@details Returns the number of species allowed on each site of the provided structure.

@param structure an atomic configuration
@param latticeSites a list of sites

@returns the number of allowed species for each site
**/
std::vector<int> ClusterSpace::getNumberOfAllowedSpeciesBySite(const Structure &structure, const std::vector<LatticeSite> &latticeSites) const
{
    std::vector<int> numberOfAllowedSpecies;
    numberOfAllowedSpecies.reserve(latticeSites.size());
    for (const auto &latsite : latticeSites)
    {
        numberOfAllowedSpecies.push_back(structure.getNumberOfAllowedSpeciesBySite(latsite.index()));
    }
    return numberOfAllowedSpecies;
}

/// Collect information about the cluster space.
/// @todo document this function
void ClusterSpace::collectClusterSpaceInfo()
{
    _clusterSpaceInfo.clear();
    std::vector<int> emptyVec = {0};
    _clusterSpaceInfo.push_back(make_pair(-1, emptyVec));
    _multiComponentVectors.resize(_orbitList.size());
    _elementPermutations.resize(_orbitList.size());
    for (int i = 0; i < _orbitList.size(); i++)
    {
        std::vector<std::vector<int>> permutedMCVector;
        auto numberOfAllowedSpecies = getNumberOfAllowedSpeciesBySite(_primitiveStructure, _orbitList.getOrbit(i).getRepresentativeSites());

        auto multiComponentVectors = _orbitList.getOrbit(i).getMultiComponentVectors(numberOfAllowedSpecies);
        if (std::none_of(numberOfAllowedSpecies.begin(), numberOfAllowedSpecies.end(), [](int n) { return n < 2; }))
        {
            auto elementPermutations = getMultiComponentVectorPermutations(multiComponentVectors, i);
            _elementPermutations[i] = elementPermutations;
            _multiComponentVectors[i] = multiComponentVectors;
        }        

        for (const auto &multiComponentVector : multiComponentVectors)
        {
            _clusterSpaceInfo.push_back(make_pair(i, multiComponentVector));
        }
    }
    _isClusterSpaceInitialized = true;
}

/// Returns information about the cluster space.
/// @todo document this function
std::pair<int, std::vector<int>> ClusterSpace::getClusterSpaceInfo(const unsigned int index)
{
    if (!_isClusterSpaceInitialized)
    {
        collectClusterSpaceInfo();
    }

    if (index >= _clusterSpaceInfo.size())
    {
        std::string msg = "Out of range in ClusterSpace::getClusterSpaceInfo: ";
        msg += std::to_string(index) + " >= " + std::to_string(_clusterSpaceInfo.size());
        throw std::out_of_range(msg);
    }

    return _clusterSpaceInfo[index];
}

//}