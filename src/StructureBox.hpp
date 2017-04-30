/** This is a testing class
    The purposes of this class are:
    1) Load a list of structures from a random regenerator
    2) Sort those cluster in a training and validation set
    2) Calculate the cluster vector of those structures
    3) Return the cluster vector and structures to ClusterExpansion and MasterJob
*/

#include "StructureLib.hpp" 
#include "ClusterBox.hpp"
#include "Structure.hpp"
#include "ClusterVector.hpp"

class StructureBox{

public:

/** Important issue: StructureLib is a library with the "DFT" structures that may also contain the primitive cell
    LatticeSet is a class that define an arrangement of atoms in a crystal structure */
    
    StructuresBox();

    ~StructureBox();
    
    std::vector<std::string> AtomsTypeSet()

    std::vector<double> CutOffSet();

    std::vector<double> Concentrations();

    /// Returns random generated atomic structures
    std::vector<Structure> StructureGeneratorInterface();

    /// Index all generated structures
    indexStructureSet();

    /// Partiotining structures set into training and validation subset
    partitioningStructureSet();

    /// Return desired subset of training or validation structure given the range of the vector
    std::vector<Structure> getStructureSet();

    /// Returns current number of structures used for training
    int getNumberStructTraining();

    /// Returns current number of strcutures used for validation
    int getNumberStructValidation();

    /// Add a structure to the training set of structures
    void addStructureTrainingSet();

    /// Set cluster functions according to the type of element and clusters
    setClusterFunctions(std::vector<std::string>& Atoms; ClusterBox::getClusters);

    /// Returns cluster vectors of running training structures and uploaded clustera in ClusterBox
    /// ClusterVector class is defined somewhere else
    std::vector<ClusterVector> ClusterVectors();

    /// Returns cluster vector of validation structures and uploaded clusters in ClusterBox
    std::vector<ClusterVector> ValidationVectors();

    /// Add new cluster vectors when training structures are added to existing set
    void addClusterVectors();

    /// Fetch cluster vector space in case of a new cutoff set
    void fetchClustersSet();
        
    /// Add new element to cluster vectors when new clusters are either added
    void addElementsClusterVector();

    // Copy only those cluster vectors and elements related to specific set of structures and clusters
    std::vector<ClusterVector> getClusterVectors();
}
