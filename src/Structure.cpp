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

PYBIND11_PLUGIN(example)
{
    py::module m("example", "pybind11 example plugin");

    py::class_<Structure>(m, "Structure")
        .def(py::init<const Eigen::Matrix<double, Dynamic, 3> &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix3d &,
                      const std::vector<bool> &>())
        .def("set_positions", &Structure::setPositions)
        .def("set_elements", &Structure::setElements)
        .def("get_elements", &Structure::getElements)
        .def("get_positions", &Structure::getPositions)
        .def("get_distance", &Structure::getDistance)
        .def("get_distance2", &Structure::getDistance2)
        .def("has_pbc", &Structure::has_pbc)
        .def("get_pbc", &Structure::get_pbc)
        .def("set_pbc", &Structure::set_pbc)
        .def("get_cell", &Structure::get_cell)
        .def("set_cell", &Structure::set_cell)
        .def("print_positions", &Structure::printPositions);
    return m.ptr();
}