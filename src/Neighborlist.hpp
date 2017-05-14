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
    std::pair<std::vector<int>, std::vector<Vector3d>> getNeighbors(int index)
    {
        return std::make_pair( latticeIndices[index], offsets[index]);
    }
    private:
    std::vector<std::vector<int>> latticeIndices;
    std::vector<std::vector<Vector3d>> offsets;        
    double _cutoff;
    double DISTTOL = 1e-7;
};
