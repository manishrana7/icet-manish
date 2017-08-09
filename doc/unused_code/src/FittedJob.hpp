/** This is a testing class
*   The purposes of this class ares:
*   1) Expose/pass the user's input to the others classes
*   2) Evaluate whether cluster expansion is dirty/ Control the number of iterations
*   3) Save progressively the outputs: rms, rms_training, fitted_values, runing_structures_index
*   4) Average ECIs
*   5) Once the expansion is clean provide ECIs, fitted_values, loocv, lpocv, runing_structures, runing_clusters
    as attributes of the class.
*   Note: This class only will load structures and clusters using the collected indexes for the final result.
*/

/*! A list of included headers may contain:*/

#include "ClusterExpansion.hpp"
#include "StructureBox.hpp"
#include "ClusterBox.hpp" 



/*! Particularly I do not like the name of the class */
class FittedJob{

public:


    FittedJob(); /// Default constructor

    ~FittedJob();

    /*! Important note: ClusterBox and StructureBox will be initilized here */

    /// Initialize desirable cluster list from cutoff vector
    void ClusterBox(std::vector<double> CutOffs);

    /// Returns the total number of structures to be generated for training
    int dimensionStructTraining();

    /// Returns the total number of structures to be generated for validation
    int dimensionStructValidation();

    /// Returns number of training structures used per round
    int rangeStructTraining();

    ///Returns number of validation structures used per round
    int rangeStructValidation();

    /// Initialize list of structures avaliable for training and validation
    void StructureBox();

    /// Returns the number of structures to be used for training in a single iteration
    void rangeStructTraining(const int numberStructTraining);

    /// Returns the number of structures to be used for validation in a single iteration
    void rangeStructValidation(const int numberStructValidation);

    /// Setting up cluster expansion by passing method and (optional) intial parameters 
    void ClusterExpansion(std::string&);

    bool _isClusterExpansionDirty=true;
 
    /// Returns whether fitting will be performed according either 
    /// to the cleaness of the expansion or the lack of avaliable  structures for training 
    /// or whatever else
    bool fittingEnabled();

    /// Load new range of structures for the next calculation
    void newRangeTraining(ClusterExpansion::rangeTrainingStruct);

    /// Returns averaged expansion coefficients over number of iteration
    std::vector<double> averagedECIs(ClusterExpansion::getECIs);

    /// Returns fitted properties 
    std::vector<double> getFittedProperties(); 

    /// Returns root-mean-square error for validation structures
    double rmsValidationError()

    /// Returns rms error for training structure
    double rmsTrainingError();

    /// Returns and load the number of iterations performed so far  
    int getNumberIterations();

    /// Returns LpOCV error
    double LPOCV();

    /// Returns LOOCV error
    double LOOCV();

    /// Returns set of structures used for training 
    std::vector<Structure>  StructureBox::getStructureSet(std::vector<int> rangeStructTraining);

    /// Returns set of structures used for validation
    std::vector<StructureLattice>  StructureBox::getStructureSet(std::vector<int> rangeStructValidation);

    /// Returns cluster vectors 
    std::vector<ClusterVector> StructureBox::getClusterVectors();

    /// Returns set of cluster used for the job
    std::vector<Structure> ClusterBox::getClusters(); 

}



