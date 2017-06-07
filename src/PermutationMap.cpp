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

    for (size_t i = 0; i < _translations.size(); i++)
    {
        std::vector<Vector3d> permutationPositionI(fractionalPositions.rows()); //do reserve instead

        for (size_t j = 0; j < fractionalPositions.size(); j++)
        {
            Vector3d permutatedPos = _translations[i] +  fractionalPositions.row(j) * _rotations[i];
            permutationPositionI[j] = permutatedPos;
        }
    }


}