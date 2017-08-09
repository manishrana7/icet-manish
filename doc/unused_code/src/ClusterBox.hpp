/** This is a testing class
    The purpose of this class are:
    1) Load a list of clusters from  Cluster Library
    2) Sort those clusters given its length
    3) Fetch difference between clusters set
*/

#include "StructureLattice.hpp"

class ClusterBox{
    
public:

    ClusterBox();

    ~ClusterBox();

    std::vector<StructureLattice> getClusters(std::vector<double>& CutOffSet);

    void indexClusters();

    int getNumberClusters();

    void addClusters();

    void fetchingClusters();

    void sortClusters();
}
