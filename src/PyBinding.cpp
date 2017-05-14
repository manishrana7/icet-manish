#include "Structure.hpp"
#include "Neighborlist.hpp"
#include "ManybodyNeighborlist.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>


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
        .def("print_positions", &Structure::printPositions)
        ;

   py::class_<Neighborlist>(m, "Neighborlist")
        .def(py::init<const double>())
        .def("build", &Neighborlist::build)
        .def("get_neighbors",&Neighborlist::getNeighbors)
        ;
           
   py::class_<ManybodyNeighborlist>(m, "ManybodyNeighborlist")
     .def(py::init<const std::vector<double>>())
     .def("calc_intersection", &ManybodyNeighborlist::getIntersection)
     ;
    return m.ptr();
}