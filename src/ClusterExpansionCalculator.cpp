#include "ClusterExpansionCalculator.hpp"

ClusterExpansionCalculator::ClusterExpansionCalculator(const ClusterSpace &clusterSpace, const Structure &structure)
{
    _clusterSpace = clusterSpace;
    _superCell = structure;
    _basisAtomOrbitList.clear();
    // Set up stuff.
    _theLog = LocalOrbitListGenerator(clusterSpace.getOrbitList(), _superCell);

    OrbitList fullOrbitList = _theLog.generateFullOrbitList();
    Vector3d zeroOffset = {0.0, 0.0, 0.0};    

    for (int i = 0; i < _clusterSpace.getPrimitiveStructure().size(); i++)
    {
        OrbitList basisOrbitList = OrbitList();
        basisOrbitList.setPrimitiveStructure(_clusterSpace.getPrimitiveStructure());

        LatticeSite primitiveSite_i = LatticeSite(i, zeroOffset);
        LatticeSite superCellEquivalent = _theLog.getPrimToSupercellMap()[primitiveSite_i];

        for (const auto orbit : fullOrbitList.getOrbitList())
        {
            Orbit basisOrbit = Orbit(orbit.getRepresentativeCluster());

            for (const auto latticeSites : orbit.getEquivalentSites())
            {
                if (std::find(latticeSites.begin(), latticeSites.end(), superCellEquivalent) != latticeSites.end())
                {
                    basisOrbit.addEquivalentSites(latticeSites);
                }
            }
            basisOrbitList.addOrbit(basisOrbit);
        }
        _basisAtomOrbitList.push_back(basisOrbitList);
    }
}


std::vector<double> ClusterExpansionCalculator::getLocalClusterVector(const Structure &structure, const int index) const
{
    bool orderIntact = true; // count the clusters in the orbit with the same orientation as the prototype cluster
    
    ClusterCounts clusterCounts = ClusterCounts();
    
    Vector3d localPosition = structure.getPositions().row(index);
    LatticeSite localSite = _clusterSpace.getOrbitList().getPrimitiveStructure().findLatticeSiteByPosition(localPosition);
    int basisIndex = localSite.index();

    clusterCounts.countOrbitList(structure, _basisAtomOrbitList[basisIndex], orderIntact);

  


    const auto clusterMap = clusterCounts.getClusterCounts();
    std::vector<double> clusterVector;
    clusterVector.push_back(1.0 / structure.size());
    // Finally begin occupying the cluster vector
    for (size_t i = 0; i < _clusterSpace.getOrbitList().size(); i++)
    {
        auto repCluster = _clusterSpace.getOrbitList().getOrbit(i).getRepresentativeCluster();
        std::vector<int> allowedOccupations;
        try { 
                allowedOccupations = _clusterSpace.getAllowedOccupations(_clusterSpace.getPrimitiveStructure(), _clusterSpace.getOrbitList().getOrbit(i).getRepresentativeSites());
            }
        catch (const std::exception& e)
        { 
            throw std::runtime_error("Failed getting allowed occupations in genereteClusterVector"); 
        }

        // Skip rest if any sites aren't active sites (i.e. allowed occupation < 2)
        if (std::any_of(allowedOccupations.begin(), allowedOccupations.end(),[](int allowedOccupation ){ return allowedOccupation < 2; }))
        {
            continue;
        }
        
        auto mcVectors = _clusterSpace.getOrbitList().getOrbit(i).getMCVectors(allowedOccupations);
        auto allowedPermutationsSet = _clusterSpace.getOrbitList().getOrbit(i).getAllowedSitesPermutations();
        auto elementPermutations = _clusterSpace.getMCVectorPermutations(mcVectors, i);
        repCluster.setClusterTag(i);
        int currentMCVectorIndex = 0;
        for (const auto &mcVector : mcVectors)
        {
            double clusterVectorElement = 0;
            int multiplicity = 0;

            for (const auto &elementsCountPair : clusterMap.at(repCluster))
            {

                // TODO check if allowedOccupations should be permuted as well.
                for (const auto &perm : elementPermutations[currentMCVectorIndex])
                {
                    auto permutedMCVector = icet::getPermutedVector(mcVector, perm);
                    auto permutedAllowedOccupations = icet::getPermutedVector(allowedOccupations, perm);
                    clusterVectorElement += _clusterSpace.getClusterProduct(permutedMCVector, permutedAllowedOccupations, elementsCountPair.first) * elementsCountPair.second;
                    multiplicity += elementsCountPair.second;
                }
            }
            clusterVectorElement /= ((double) multiplicity * structure.size() );
            clusterVector.push_back(clusterVectorElement * mcVector.size());

            currentMCVectorIndex++;
        }
    }
    return clusterVector;

} 
