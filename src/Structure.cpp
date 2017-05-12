#include "Structure.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
namespace py = pybind11;

Structure::Structure()
{
    //empty constructure..
}

int add(int i, int j) {
    return i + j;
}


PYBIND11_PLUGIN(example)
{
    py::module m("example", "pybind11 example plugin");

    py::class_<Structure> structure(m, "Structure");
    structure        
        .def(py::init<>())
        .def("printHello", &Structure::printHello);

    m.def("add", &add, "A function which adds two numbers");    
 return m.ptr();
}