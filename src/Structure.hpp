#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
using namespace Eigen;

namespace py = pybind11;

class Structure
{
  public:
    Structure(const Eigen::Matrix<float, Dynamic, 3> &pos,
             const std::vector<std::string> elements);

    void setPositions(const Eigen::Matrix<float, Dynamic, 3> &);
    void setElements(const std::vector<std::string> &);
    Eigen::Matrix<float, Dynamic, 3> &getPositions();
    double getDistance(const int, const int, const bool) const;
    void printPositions();

  private:
    Eigen::Matrix<float, Dynamic, 3> positions;
    std::vector<std::string> _elements;
};
