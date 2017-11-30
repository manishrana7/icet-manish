#include "Symmetry.hpp"

namespace icet
{

/// Return the transformed position `position` using the input translation and rotation
Eigen::Vector3d transformPosition2(const Eigen::Vector3d &position, const Eigen::Vector3d &translation, const Eigen::Matrix3d &rotation)
{
    Eigen::Vector3d transformedPosition = position; //translation.transpose() + position.transpose() * rotation.transpose();
    return transformedPosition;
}

bool next_cartesian_product(const std::vector<std::vector<int>> &items, std::vector<int> &currentProduct)
{
    auto n = items.size();
    if (n != currentProduct.size())
    {
        throw std::runtime_error("ERROR: items and currentProduct are different sizes in function next_cartesian_product");
    }

    for (size_t i = 0; i < n; ++i)
    {
        if (++currentProduct[i] == items[i].size())
        {
            currentProduct[i] = 0;
        }
        else
        {
            return true;
        }
    }
    return false;
}

Eigen::Matrix3i getUnitcellPermutation(const Eigen::Matrix3d &inputCell, const Eigen::Matrix3d &referenceCell, double tolerance_cell)
{
    // obtain the (in general non-integer) transformation matrix
    // connecting the input configuration to the reference structure
    // L = L_p.P --> P = L_p^-1.L
    Eigen::Matrix3d permutationMatrix = inputCell * referenceCell.inverse();

    Eigen::Matrix3d roundedPermutationMatrix;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            roundedPermutationMatrix(i, j) = roundedPermutationMatrix(i, j);
        }
    }

    Eigen::Matrix3i integerPermutationMatrix = roundedPermutationMatrix.cast<int>();
    Eigen::Matrix3d deviationOfCells = permutationMatrix - integerPermutationMatrix.cast<double>();

    if (deviationOfCells.sum() / 9 > tolerance_cell)
    {
        std::string errorMessage = "Failed to map input cell to reference cell\n";
        errorMessage += "tolerance_cell: " + std::to_string(tolerance_cell) + " \n";
        errorMessage += "Deviation: " + std::to_string(deviationOfCells.sum() / 9) + " \n";
        errorMessage += "You can try to raise tolerance_cell\n";
        throw std::runtime_error(errorMessage);
    }

    return integerPermutationMatrix;
}
}
