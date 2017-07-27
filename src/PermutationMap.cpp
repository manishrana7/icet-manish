#include "PermutationMap.hpp"

/**
Will create all the permutated positions from
these positions and the current rotational/translational symmetries.

@TODO: Think about pruning positions that fall outside cell if pbc is false
@TODO: Relate positions to indices
@TODO: Think about possibility to sort permutationmap (both row-wise and col-wise?)
@TODO: Think about possibility to only add permutations that are bigger/smaller with motivation
       of removing duplicates. 


*/

void PermutationMap::build(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &fractionalPositions)
{
    _permutatedPositions.clear();
    _permutatedPositions.resize(fractionalPositions.rows());
    for (size_t j = 0; j < fractionalPositions.rows(); j++) //row
    {
        for (size_t i = 0; i < _translations.size(); i++) //column
        {
            Vector3d permutatedPos = _translations[i].transpose() + fractionalPositions.row(j) * _rotations[i]; // transpose frac pos?
            roundVector3d(permutatedPos);
            _permutatedPositions[j].push_back(permutatedPos);
        }
    }
}
