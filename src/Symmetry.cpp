#include "Symmetry.hpp"

namespace icet
{


/// Return the transformed position `position` using the input translation and rotation
Eigen::Vector3d transformPosition2(const Eigen::Vector3d &position, const Eigen::Vector3d &translation, const Eigen::Matrix3d &rotation)
{
    Eigen::Vector3d transformedPosition = position; //translation.transpose() + position.transpose() * rotation.transpose();
    return transformedPosition;
}
}