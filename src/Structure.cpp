#include "Structure.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
using namespace Eigen;

namespace py = pybind11;

Structure::Structure()
{
    //empty constructure..
}

void Structure::setPositions(Eigen::Matrix<float, Dynamic, 3>  &m)
{
  position = m;
}
 Eigen::Matrix<float, Dynamic, 3> & Structure::getPositions()
{
  return  position;
}
void Structure::printPositions()
{
  std::cout<<position<<std::endl;
}




PYBIND11_PLUGIN(example)
{
    py::module m("example", "pybind11 example plugin");

    py::class_<Structure> (m, "Structure")    
        .def(py::init<>())
        .def("set_positions",&Structure::setPositions)
        .def("get_positions",&Structure::getPositions)
        .def("print_position",&Structure::printPositions)
         ;   
 return m.ptr();
}