#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
using namespace Eigen;

namespace py = pybind11;

class Structure {
public:
    Structure();
    void setPositions(Eigen::Matrix<float, Dynamic, 3> &);
  Eigen::Matrix<float, Dynamic, 3>& getPositions();
  void printPositions();  
 private:
  Eigen::Matrix<float, Dynamic, 3> position;  
};


