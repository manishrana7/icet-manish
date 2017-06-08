#include "PermutationMap.hpp"

/**
Will create all the permutated positions from
these positions and the current rotational/translational symmetries.

@TODO: Think about pruning positions that fall outside cell if pbc is false
@TODO: Relate positions to indices

*/

void PermutationMap::build(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &fractionalPositions)
{
    //std::vector<std::vector<Vector3d>> _permutatedPositons;
    std::cout << "trans size " << _translations.size() << " rot size " << _rotations.size() << std::endl;

    for (size_t i = 0; i < _translations.size(); i++)
    {
        std::vector<Vector3d> permutationPositionI; //do reserve instead
        permutationPositionI.reserve(fractionalPositions.rows());
        for (size_t j = 0; j < fractionalPositions.rows(); j++)
        {
            Vector3d permutatedPos = _translations[i].transpose() + fractionalPositions.row(j).transpose() * _rotations[i]; // transpose frac pos?
            permutationPositionI.push_back(permutatedPos);
        }
        _permutatedPositions.push_back(permutationPositionI);
    }
}