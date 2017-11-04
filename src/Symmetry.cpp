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
}
