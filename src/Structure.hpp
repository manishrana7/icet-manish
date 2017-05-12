#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
namespace py = pybind11;

class Structure {
public:
    Structure();
    void printHello(int n)
    {
        for(int i=0; i<n; i++)
        {
            std::cout<<"Hello"<<std::endl;
        }
    }

};


