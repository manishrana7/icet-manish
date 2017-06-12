#include "Structure.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
using namespace Eigen;
namespace py = pybind11;

Structure::Structure(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &pos,
                     const std::vector<std::string> &elements,
                     const Eigen::Matrix3d &cell,
                     const std::vector<bool> &pbc)
{
    setPositions(pos);
    setStrElements(elements);    
    _cell = cell;
    _pbc = pbc;
    _uniqueSites.resize(elements.size());
}

/**

Returns the distance between two indices with MIC = false
@todo:  add  mic functionality

*/

double Structure::getDistance(const int index1, const int index2) const
{
    

    if (index1 >= _positions.rows() or index2 >= _positions.rows())
    {
        throw std::out_of_range("Error: Tried accessing position at out of bound index. Structure::getDistance");
    }

    const Eigen::Vector3d pos1 = _positions.row(index1);
    const Eigen::Vector3d pos2 = _positions.row(index2);

    return (pos1 - pos2).norm();
}
