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
    Structure(const Eigen::Matrix<double, Dynamic, 3> &,
              const std::vector<std::string> &,
              const Eigen::Matrix<double, 3, 3> &,
              const std::vector<bool> &);

    double getDistance(const int, const int, const bool) const;
    void printPositions();

    // Getters - Setters
    void setPositions(const Eigen::Matrix<double, Dynamic, 3> &positions)
    {
        _positions = positions;
    }
 

    Eigen::Matrix<double, Dynamic, 3> &getPositions()
    {
        return _positions;
    }

   void setElements(const std::vector<std::string> &elements)
    {
        _elements = elements;
    }

    std::vector<std::string> getElements() const
    {
        return _elements;
    }
    bool has_pbc(const int k) const
    {
        return _pbc[k];
    }
    std::vector<bool> get_pbc() const
    {
        return _pbc;
    }
    void set_pbc(const std::vector<bool> pbc)
    {
        _pbc = pbc;
    }

    void set_cell( const Eigen::Matrix<double, 3, 3> &cell)
    {
        _cell = cell;
    }

   Eigen::Matrix<double, 3, 3> get_cell() const
    {
        return _cell;
    }


  private:
    Eigen::Matrix<double, Dynamic, 3> _positions;
    Eigen::Matrix<double, 3, 3> _cell;
    std::vector<std::string> _elements;
    std::vector<bool> _pbc;
};
