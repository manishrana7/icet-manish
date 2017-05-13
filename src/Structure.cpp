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

Structure::Structure(const Eigen::Matrix<double, Dynamic, 3>  &pos,
const std::vector<std::string> elements)
{
    setPositions(pos);
    setElements(elements);
}

void Structure::setPositions(const Eigen::Matrix<double, Dynamic, 3>  &m)
{
  positions = m;
}

void Structure::setElements(const std::vector<std::string> &elements)
{
  _elements = elements;
}

double Structure::getDistance(const int index1, const int index2, const bool mic) const
{
    const Eigen::Vector3d pos1 = positions.row(index1);
    const Eigen::Vector3d pos2 = positions.row(index2);

    Eigen::Vector3d diff = pos1 - pos2;
    if(mic)
    {

    }

   const  auto distance = diff.norm();
    return distance;
    

}

 Eigen::Matrix<double, Dynamic, 3> & Structure::getPositions() 
{
  return  positions;
}
void Structure::printPositions() 
{
  std::cout<<positions<<std::endl;
}




PYBIND11_PLUGIN(example)
{
    py::module m("example", "pybind11 example plugin");

    py::class_<Structure> (m, "Structure")    
        .def(py::init<const Eigen::Matrix<double, Dynamic, 3>  &,
    const std::vector<std::string> >()) 
        .def("set_positions",&Structure::setPositions)
        .def("set_elements", &Structure::setElements)
        .def("get_positions",&Structure::getPositions)
        .def("get_distance",&Structure::getDistance)
        .def("print_positions",&Structure::printPositions)
         ;   
 return m.ptr();
}