#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <string>
#include "Structure.hpp"

using namespace Eigen;

/**
class Neighborlist

Builds a (pair) neighborlist for each lattice site

*/

class Neighborlist
{
    public:
    Neighborlist(const double);
    void build(const Structure &);
    void update(const Structure &);
    // std::pair<std::vector<int>, std::vector<Vector3d>> getNeighbors(int index)
    // {
    //     return std::make_pair( latticeIndices[index], offsets[index]);
    // }
    std::vector<std::pair<int,Vector3d>> getNeighbors(int index) const
    {
        return _neighbors[index];
    }
    ///Check if index1 and index2 are neighbors
    bool isNeighbor(const int index1, const int index2) const
    {
        for(const auto &nbr1 : _neighbors[index1])
        {
            if( nbr1.first == index2 )
            {
                return true;
            }
        }
    }

     size_t size() const
    {
        return _neighbors.size();
    }
    private:
    std::vector<std::vector<int>> latticeIndices;
    std::vector<std::vector<Vector3d>> offsets;
    std::vector<std::vector<std::pair<int,Vector3d>>> _neighbors;        
    double _cutoff;
    double DISTTOL = 1e-7;
};
