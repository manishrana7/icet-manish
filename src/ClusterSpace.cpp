#define _USE_MATH_DEFINES
#include <cmath>

#include "ClusterSpace.hpp"

/**
@details This constructor initializes a ClusterSpace object.
@param chemicalSymbols vector of allowed chemical symbol for each site
@param orbitList list of orbits for the primitive structure
@param positionTolerance
    tolerance applied when comparing positions in Cartesian coordinates
@param fractionalPositionTolerance
    olerance applied when comparing positions in fractional coordinates
*/
ClusterSpace::ClusterSpace(std::shared_ptr<OrbitList> orbitList,
                           const double positionTolerance,
                           const double fractionalPositionTolerance)
    : _primitiveOrbitList(orbitList)
{
    // Set up a map between atomic numbers and the internal species enumeration scheme.
    for (const auto &atomicNumbers : _primitiveOrbitList->structure().allowedAtomicNumbers())
    {
        std::unordered_map<int, int> speciesMap;
        std::vector<int> atomicNumbersCopy = atomicNumbers;
        sort(atomicNumbersCopy.begin(), atomicNumbersCopy.end());
        for (size_t i = 0; i < atomicNumbersCopy.size(); i++)
        {
            speciesMap[atomicNumbersCopy[i]] = i;
        }
        _speciesMaps.push_back(speciesMap);
    }
}

/**
@details This function calculates and then returns the cluster vector for the
input structure in the cluster space.

The first element in the cluster vector will always be one (1) corresponding to
the zerolet. The remaining elements of the cluster vector represent averages
over orbits (symmetry equivalent clusters) of increasing order and size.

@param structure
    Input configuration
@param fractionalPositionTolerance
    Tolerance applied when comparing positions in fractional coordinates
**/

std::vector<double> ClusterSpace::getClusterVector(const Structure &structure,
                                                   const double fractionalPositionTolerance) const
{
    // Construct orbit list for this structure.
    std::shared_ptr<Structure> supercell = std::make_shared<Structure>(structure);
    LocalOrbitListGenerator localOrbitListGenerator = LocalOrbitListGenerator(*_primitiveOrbitList,
                                                                              supercell,
                                                                              fractionalPositionTolerance);
    auto currentOrbitList = localOrbitListGenerator.getFullOrbitList();

    // Check that the number of unique offsets equals the number of unit cells in the supercell.
    if (localOrbitListGenerator.getNumberOfUniqueOffsets() !=
        structure.size() / _primitiveOrbitList->structure().size())
    {
        std::ostringstream msg;
        msg << "The number of unique offsets does not match the number of "
            << "primitive units in the input structure(ClusterSpace::getClusterVector) "
            << std::endl
            << localOrbitListGenerator.getNumberOfUniqueOffsets() << " != "
            << structure.size() / _primitiveOrbitList->structure().size();
        throw std::runtime_error(msg.str());
    }
    return getClusterVectorFromOrbitList(currentOrbitList, supercell);
}

/**
@details Evaluates the cluster function using the specified parameters.

The cluster functions (also "orthogonal point functions") are defined as

.. math::

   \\Theta_{n}(\\sigma_p) = \\begin{cases}
      1                                     &\\quad \\text{if}~n=0 \\
      -\\cos\\left(\\pi(n+1)\\sigma_p/M\\right) &\\quad \\text{if n is odd} \\
      -\\sin\\left(\\pi n   \\sigma_p/M\\right) &\\quad \\text{if n is even}
    \\end{cases}

@param numberOfAllowedSpecies number of allowed species on the site in question
@param clusterFunction index of cluster function
@param species index of species

@returns the value of the cluster function
*/
double ClusterSpace::evaluateClusterFunction(const int numberOfAllowedSpecies,
                                             const int clusterFunction,
                                             const int species) const
{
    if (((clusterFunction + 2) % 2) == 0)
    {
        return -cos(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) *
                    (double)species / ((double)numberOfAllowedSpecies));
    }
    else
    {
        return -sin(2.0 * M_PI * (double)((int)(clusterFunction + 2) / 2) *
                    (double)species / ((double)numberOfAllowedSpecies));
    }
}

/**
@details Evaluates the full cluster product of the entire cluster.

@param multicomponentVector
    Multicomponent vector, each element of the vector gives the index of a
    cluster function
@param numberOfAllowedSpecies
    Number of species allowed on the sites in this cluster.
@param species
    Species that occupy (decorate) the cluster identified by atomic number
@param indices
    Representative lattice indices of the cluster being computed
@param permutation
    Describes the desired order of the points in the cluster; for example,
    if the current cluster is a triplet cluster and permutation = {1, 0, 2},
    then multiComponentVector, numberOfAllowedSpecies, and indices should
    be permuted such that their first two elements change places

@returns the cluster product
**/
double ClusterSpace::evaluateClusterProduct(const std::vector<int> &multicomponentVector,
                                            const std::vector<int> &numberOfAllowedSpecies,
                                            const std::vector<int> &species,
                                            const std::vector<int> &indices,
                                            const std::vector<int> &permutation) const
{
    double clusterProduct = 1;

    for (size_t i = 0; i < species.size(); i++)
    {
        int index = permutation[i];
        clusterProduct *= evaluateClusterFunction(numberOfAllowedSpecies[index],
                                                  multicomponentVector[index],
                                                  _speciesMaps[indices[index]].at(species[i]));
    }
    return clusterProduct;
}

/**
@brief Calculates the size of the cluster space, defined as the length of a cluster vector.
*/
size_t ClusterSpace::size() const
{
    int size = 1; // 1 for the zerolet
    for (size_t orbitIndex = 0; orbitIndex < _primitiveOrbitList->size(); orbitIndex++)
    {
        const Orbit &orbit = _primitiveOrbitList->getOrbit(orbitIndex);
        size += orbit._clusterVectorElements.size();
    }
    return size;
}

/*
@details
    Occupy cluster vector based on a supercell and a corresponding orbit list.

@param orbitList
    An orbit list to be used for counting. Can be either a full orbit list
    or a "local orbit list"
@param supercell
    Defines the occupations of the structure whose cluster vector should be
    computed.
@param firstElement
    First element of the cluster vector (default: 1.0).
@param flipIndex
    If a local cluster vector should be calculated this argument is used to
    specify the index of the site whose local cluster vector should be
    computed. If calculating a change in cluster vector, this is the site
    whose occupation has changed. If -1 (default), the total cluster vector
    will be calculated.
@param newOccupation
    New atomic number on the site with index flipIndex. If this argument is
    not -1, a change in cluster vector will be calculated.
*/
const std::vector<double> ClusterSpace::getClusterVectorFromOrbitList(const OrbitList &orbitList,
                                                                      const std::shared_ptr<Structure> supercell,
                                                                      const double firstElement,
                                                                      const int flipIndex,
                                                                      const int newOccupation) const
{
    if (newOccupation >= 0 && flipIndex == -1)
    {
        std::ostringstream msg;
        msg << "flipIndex needs to be specified (larger than -1) if newOccupation "
            << "is specified (ClusterSpace::getClusterVectorFromOrbitList)";
        throw std::runtime_error(msg.str());
    }
    std::vector<double> clusterVector;
    clusterVector.push_back(firstElement);
    if (_primitiveOrbitList->size() != orbitList.size())
    {
        std::ostringstream msg;
        msg << "Orbit lists do not match (ClusterSpace::getClusterVectorFromOrbitList)."
            << std::endl
            << orbitList.size() << " >= " << _primitiveOrbitList->size() << std::endl;
        throw std::runtime_error(msg.str());
    }

    // Define some variable before loop for performance reasons
    std::map<std::vector<int>, double> clusterCounts;
    std::vector<int> allowedOccupations;
    std::vector<int> indicesOfRepresentativeSites;

    // Start to occupy the cluster vector orbit by orbit
    for (size_t orbitIndex = 0; orbitIndex < _primitiveOrbitList->size(); orbitIndex++)
    {
        const Orbit &supercellOrbit = orbitList.getOrbit(orbitIndex);
        const Orbit &primitiveOrbit = _primitiveOrbitList->getOrbit(orbitIndex);

        // Count clusters
        if (newOccupation > -1)
        {
            clusterCounts = supercellOrbit.getClusterCountChanges(supercell, flipIndex, newOccupation);
        }
        else
        {
            clusterCounts = supercellOrbit.getClusterCounts(supercell, flipIndex);
        }

        // Extract allowed occupations
        allowedOccupations = primitiveOrbit.representativeCluster().getNumberOfAllowedSpeciesPerSite();

        // Skip the rest if any of the sites are inactive (i.e. allowed occupation < 2)
        // @todo Let orbit handle this
        if (std::any_of(allowedOccupations.begin(), allowedOccupations.end(), [](int allowedOccupation)
                        { return allowedOccupation < 2; }))
        {
            continue;
        }

        indicesOfRepresentativeSites.clear();
        for (const LatticeSite &site : primitiveOrbit.representativeCluster().latticeSites())
        {
            indicesOfRepresentativeSites.push_back(site.index());
        }

        // Loop over all multicomponent vectors for this orbit.
        // These are vectors of integers (where the integer represents a cluster function index).
        //
        // Example 1: For a binary alloy we obtain [0, 0] and [0, 0, 0]
        // for pair and triplet terms, respectively.
        //
        // Example 2: For a ternary alloy we obtain [0, 0], [0, 1], [1, 1] for pairs
        // and similarly for triplets.
        //
        // Depending on the symmetry of the cluster one might also obtain [1, 0]
        // (e.g., in a clathrate or for some clusters on a HCP lattice).
        for (auto &cvElement : primitiveOrbit.clusterVectorElements())
        {
            double clusterVectorElement = 0;

            /// Loop over all the counts for this orbit
            for (const auto &clusterCount : clusterCounts)
            {
                /// Loop over all equivalent permutations for this orbit and mc vector
                for (const auto &perm : cvElement.sitePermutations)
                {
                    clusterVectorElement += clusterCount.second *
                                            evaluateClusterProduct(cvElement.multicomponentVector,
                                                                   allowedOccupations,
                                                                   clusterCount.first,
                                                                   indicesOfRepresentativeSites,
                                                                   perm);
                }
            }

            // Usually we could have counted multiplicity by simply adding the number of
            // clusters in the orbit (elementsCountPair.second), but in the case of
            // local cluster vectors or changes in cluster vectors, we have only counted
            // a subset of the clusters. We therefore use the pre-computed multiplicity.
            clusterVector.push_back(clusterVectorElement / (double)cvElement.multiplicity *
                                    (double)_primitiveOrbitList->structure().size() / (double)supercell->size());
        }
    }
    return clusterVector;
}
