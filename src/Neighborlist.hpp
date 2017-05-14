#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
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
    private:
    std::vector<std::vector<int>> latticeIndices;
    std::vector<std::vector<Vector3d>> offsets;        
    double _cutoff;
    double DISTTOL = 1e-7;
};
