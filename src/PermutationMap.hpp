#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <pybind11/stl.h>
#include <map>
#include <unordered_map>
#include <algorithm>

using namespace Eigen;

namespace py = pybind11;

class PermutationMap
{
  public:
    PermutationMap(const std::vector<Vector3d> &translations,
                   const std::vector<Matrix3d> &rotations)
    {
        _translations = translations;
        _rotations = rotations;
    }

    void build(const Eigen::Matrix<double,Dynamic,3,RowMajor> &positions);

  private:
    std::vector<Vector3d> _translations;
    std::vector<Matrix3d> _rotations;
    std::vector<std::vector<Vector3d>> _permutatedPositons;
    
    };