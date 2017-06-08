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

    void build(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &positions);

    std::vector<std::vector<Vector3d>> getPermutatedPositions()
    {
        return _permutatedPositions;
    }


    /**

     Returns indices for unique positions,
     same indices share the same indice

    */
    std::vector<std::vector<int>> getIndicedPermutatedPositions()
    {   
        std::vector<Vector3d> uniquePositions;
        std::vector<std::vector<int>> indicePositions;
        for(const auto &posRow : _permutatedPositions)
        {
            std::vector<int> indiceRow(posRow.size());
            int indiceCount = 0;
            for(const Vector3d &pos : posRow)
            {
                const auto find = std::find(uniquePositions.begin(), uniquePositions.end(), pos);
                if(find == uniquePositions.end())
                {
                    uniquePositions.push_back(pos);
                    indiceRow[indiceCount++] = uniquePositions.size()-1;
                }
                else
                {
                    indiceRow[indiceCount++] = std::distance(uniquePositions.begin(), find);
                }
            }
            indicePositions.push_back(indiceRow);
        }

        return indicePositions;
    }

  private:
    std::vector<Vector3d> _translations;
    std::vector<Matrix3d> _rotations;
    std::vector<std::vector<Vector3d>> _permutatedPositions;
};