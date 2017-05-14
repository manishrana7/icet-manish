#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include "Vector3dCompare.hpp"
#include <vector>

/**
Design approach:
    input pair neighbors and calculate higher order neighbors
    using set intersection.
*/

class ManybodyNeighborlist
{
    public:
    ManybodyNeighborlist(const std::vector<double> cutoffs)
    {
        _cutoffs = cutoffs;
    }
    void build();

std::vector<std::pair<int,Vector3d>> getIntersection(const std::vector<std::pair<int,Vector3d>> &Ni,const std::vector<std::pair<int,Vector3d>> &Nj)
{
  std::vector<std::pair<int,Vector3d>> N_intersection;
  std::set_intersection(Ni.begin(), Ni.end(),
                            Nj.begin(), Nj.end(),
                            std::back_inserter(N_intersection),
                            NeighborPairCompare());
  return N_intersection;
}



    private:
    std::vector<double> _cutoffs;

};