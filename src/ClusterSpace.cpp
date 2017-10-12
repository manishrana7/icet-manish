#include "ClusterSpace.hpp"

/**
Calculate a cluster vector from the structure and  the internal state of the clusterspace.

The first element in the clustervector will always be a constant (1)

The second element in the clustervector is the average of the first orbit.
The third element might be the first orbits second multicomponent vector
or it is the second orbit if it only has one multicomponent vector.

The length of the clustervector will be 1 + sum_i( orbit_i * orbit_i.multicimponentvectors.size() )

*/
std::vector<double> ClusterSpace::generateClustervector(const Structure &structure) const
{

    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    LocalOrbitlistGenerator localOrbitListGenerator = LocalOrbitlistGenerator(_primitive_orbitlist, structure);
    size_t uniqueOffsets = localOrbitListGenerator.getUniqueOffsetsCount();

    ClusterCounts clusterCounts = ClusterCounts();
    for (int i = 0; i < uniqueOffsets; i++)
    {
        const auto local_orbitlist = localOrbitListGenerator.generateLocalOrbitlist(i);
        clusterCounts.countOrbitlist(structure, local_orbitlist, orderIntact);
    }

    const std::unordered_map<Cluster, std::map<std::vector<int>, int>> clusterMap = clusterCounts.getClusterCounts();
    std::vector<double> clusterVector;
    clusterVector.push_back(1);
    //Finally begin occupying the clustervector
    for (size_t i = 0; i < _primitive_orbitlist.size(); i++)
    {
        auto repCluster = _primitive_orbitlist.getOrbit(i).getRepresentativeCluster();
        repCluster.setClusterTag(i);

        double clusterVectorElement = 0;
        int multiplicity = 0;
        for (const auto &elementsCountPair : clusterMap.at(repCluster))
        {
            clusterVectorElement += getClusterProduct(_Mi, elementsCountPair.first) * elementsCountPair.second;
            multiplicity += elementsCountPair.second;
        }
        clusterVectorElement /= ((double)multiplicity);
        clusterVector.push_back(clusterVectorElement);
    }
    return clusterVector;
}

/**
This is the default clusterfunction
*/
double ClusterSpace::defaultClusterFunction(const int Mi, const int element) const
{
    // std::cout<<Mi<< " "<< element<< " = "<<(element % 2) * 2 - 1<<std::endl;

    return (element % 2) * 2 - 1;
}

///Return the full cluster product of entire cluster (elements vector). Assuming all sites have same Mi
double ClusterSpace::getClusterProduct(const int Mi, const std::vector<int> &elements) const
{
    double clusterProduct = 1;
    for (const auto &element : elements)
    {
        // std::cout<<clusterProduct<< " "<< defaultClusterFunction(Mi, element)<< std::endl;
        clusterProduct *= defaultClusterFunction(Mi, element);
    }
    return clusterProduct;
}