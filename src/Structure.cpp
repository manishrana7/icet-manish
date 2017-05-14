#include "Structure.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace py = pybind11;

Structure::Structure(const Eigen::Matrix<double, Dynamic, 3> &pos,
                     const std::vector<std::string> &elements,
                     const Eigen::Matrix3d &cell,
                     const std::vector<bool> &pbc)
{
    setPositions(pos);
    setElements(elements);
    _cell = cell;
    _pbc = pbc;
}

double Structure::getDistance(const int index1, const int index2, const bool mic) const
{
    const Eigen::Vector3d pos1 = _positions.row(index1);
    const Eigen::Vector3d pos2 = _positions.row(index2);

    Eigen::Vector3d diff = pos1 - pos2;
    if (mic)
    {
        // insert mic finder function
    }

    const auto distance = diff.norm();
    return distance;
}

void Structure::printPositions()
{
    std::cout << _positions << std::endl;
}

