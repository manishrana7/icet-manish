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

    Eigen::Matrix3d roundedPermutationMatrix = permutationMatrix;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            roundedPermutationMatrix(i, j) = round(roundedPermutationMatrix(i, j));
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

std::vector<Eigen::Matrix3i> getUnitcellSubPermutations(Eigen::Matrix3i &permutation_matrix)
{
    // 1 get vector form of matrix
    Eigen::Map<Eigen::RowVectorXi> v1(permutation_matrix.data(), permutation_matrix.size());

    // Create the cartesian items that will be used to generate sub permutations
    std::vector<std::vector<int>> cartesianItems(9, std::vector<int>(1));
    for (int i = 0; i < 9; i++)
    { //In case an element is negative,
        int sign_of_vi = 1;
        if (v1(i) != 0)
        {
            sign_of_vi = v1(i) / abs(v1(i));
        }
        for (int j = 1; j < abs(v1(i)); j++)
        {
            // Multiply with sign as to add cartesian items from -inf to zero
            cartesianItems[i].push_back(j * sign_of_vi);
        }
        std::sort(cartesianItems[i].begin(), cartesianItems[i].end());
    }

    // 2 Create vector representation of subpermutation matrices
    std::vector<int> currentMatrix(9);
    for (int i = 0; i < 9; i++)
    {
        currentMatrix[i] = cartesianItems[i][0];
    }


    // 3 Generate the sub permutation matrices with cartesian product
    std::vector<Eigen::Matrix3i> subPermutationMatrices;
    do
    {
        Eigen::Matrix3i matrixSubPermutation(currentMatrix.data());
        subPermutationMatrices.push_back(matrixSubPermutation);

    } while (next_cartesian_product(cartesianItems, currentMatrix));


    return subPermutationMatrices;
}

/// Creates a vector of positions by offsetting the structure's positions by an offset
std::vector<Eigen::Vector3d> getOffsetPositions(const Structure &structure, const Eigen::Vector3d &offset)
{
    std::vector<Eigen::Vector3d> positions;
    positions.reserve(structure.size());
    for(int i=0; i<structure.size(); i++)
    {
        LatticeSite latticeSite = LatticeSite(i, offset);
        Eigen::Vector3d position = structure.getPosition(latticeSite);
        positions.push_back(position);
    }
    return positions;

}


}
