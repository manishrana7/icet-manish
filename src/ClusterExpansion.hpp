/** This is a testing class
    The purpose of this class are:
    1) Calculated coefficient for a single round of single expansion 
    2) Use a desired compressive sampling technique or use split-Bregman by default.
    3) Calculated fitted and predicted properties from current cluster vectors
*/

#include "StructureBox.hpp"


class ClusterExpansion{

public:

    ClusterExpansion();
  
    ~ClusterExpansion();

    std::string CompressiveSampling = "split-bregman";

    int numberStructTraining;

    int numberStructValidation;

    void initMu();

    void initLambda();
    
    /// Returns a vector of cluster vectors running fot the cluster expansion
    std::vector<ClusterVector> ClusterVectors StructureBox::getClusterVectors(FittedJob::rangeStructTraining);

    /// Returns a vector of validation vector running for the cross validation
    std::vector<ClusterVector> ValidationVectors StructureBox::getClusterVectors(FittedJob::rangeStructValidation);

    /// Returns vector with the value of true properties for training
    std::vector<double> getTrueProperties();

    /// Returns the expansion coefficients according to the compressive sampling method desired
    /// Those methods may be defined somewhere else
    std::vector<double> getECIs(std::vector<ClusterVector>& , std::vector<double>&, std::string);

    /// Returns a vector with the values fitted properties
    std::vector<double> getFittedProperties(std::vector<ClusterVector>&, std::vector<double>&);
    
    /// Returns a vector with the values of the predidted properties
    /// obtained from the averaged ECIs
    std::vector<double> getPredictedProperties(std::vector<ClusterVector>&, std::vector<double>&);

    /// Deallocate memory once round calculation is done
    ClusterVectors.clear();

    /// Deallocate memory once job is finished
    ValidationVector.clear();

}
