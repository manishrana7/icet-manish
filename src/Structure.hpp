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
    Structure(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &,
              const std::vector<std::string> &,
              const Eigen::Matrix3d &,
              const std::vector<bool> &);

    double getDistance(const int, const int, const bool) const;

    double getDistance2(const int ind1, const Vector3d offset1,
                        const int ind2, const Vector3d offset2) const
    {
        Vector3d pos1 = _positions.row(ind1) + offset1.transpose() * _cell;
        Vector3d pos2 = _positions.row(ind2) + offset2.transpose() * _cell;
        Vector3d d = pos1 - pos2;
        return d.norm();
    }

    void printPositions();

    // Getters - Setters
    void setPositions(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &positions)
    {
        _positions = positions;
    }

    Eigen::Matrix<double, Dynamic, 3, RowMajor> &getPositions()
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

    void set_cell(const Eigen::Matrix<double, 3, 3> &cell)
    {
        _cell = cell;
    }

    Eigen::Matrix<double, 3, 3> get_cell() const
    {
        return _cell;
    }

  private:
    Eigen::Matrix<double, Dynamic, 3, RowMajor> _positions;
    Eigen::Matrix3d _cell;
    std::vector<std::string> _elements;
    std::vector<bool> _pbc;
};
