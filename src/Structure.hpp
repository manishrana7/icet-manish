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

    double getDistance(const int, const int) const;

    /**
        Returns the distance for index1 with unitcell offset offset 1 to index2 with unit-cell offset offset2
    */
    double getDistance2(const int index1, const Vector3d offset1,
                        const int index2, const Vector3d offset2) const
    {
        if (index1 >= _positions.rows() or index2 >= _positions.rows())
        {
            throw std::out_of_range("Error: Tried accessing position at out of bound index. Structure::getDistance2");
        }

        Vector3d pos1 = _positions.row(index1) + offset1.transpose() * _cell;
        Vector3d pos2 = _positions.row(index2) + offset2.transpose() * _cell;

        return (pos1 - pos2).norm();
    }
   
    // Getters - Setters
    void setPositions(const Eigen::Matrix<double, Dynamic, 3> &positions)
    {
        _positions = positions;
    }

    Eigen::Matrix<double, Dynamic, 3, RowMajor> getPositions() const
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

    size_t size() const
    {
        if (_elements.size() != _positions.rows())
        {
            throw std::out_of_range("Error: Positions and elements do not match in size");
        }        
        return( _elements.size());
    }    

  private:
    Eigen::Matrix<double, Dynamic, 3, RowMajor> _positions;
    Eigen::Matrix3d _cell;
    std::vector<std::string> _elements;
    std::vector<bool> _pbc;
};
