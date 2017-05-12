#include "Structure.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
using namespace Eigen;

namespace py = pybind11;

Structure::Structure(const Eigen::Matrix<float, Dynamic, 3>  &pos,
const std::vector<std::string> elements)
{
    setPositions(pos);
    setElements(elements);
}

void Structure::setPositions(const Eigen::Matrix<float, Dynamic, 3>  &m)
{
  position = m;
}

void Structure::setElements(const std::vector<std::string> &elements)
{
  _elements = elements;
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
        .def(py::init<const Eigen::Matrix<float, Dynamic, 3>  &,
const std::vector<std::string> >()) 
        .def("set_positions",&Structure::setPositions)
        .def("set_elements", &Structure::setElements)
        .def("get_positions",&Structure::getPositions)
        .def("print_position",&Structure::printPositions)
         ;   
 return m.ptr();
}