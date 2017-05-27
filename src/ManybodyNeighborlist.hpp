#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include "Vector3dCompare.hpp"
#include "Neighborlist.hpp"
#include <vector>

/**
Design approach:
    input pair neighbors and calculate higher order neighbors
    using set intersection.
*/

class ManybodyNeighborlist
{
  public:
    ManybodyNeighborlist()
    {
        //empty...
    }

    std::vector<std::pair<std::vector<std::pair<int, Vector3d>>,  std::vector<std::pair<int, Vector3d>>>> build(const Neighborlist &nl, int index, int order,bool);

    void combineToHigherOrder(const Neighborlist &nl, 
                                 std::vector<std::pair<std::vector<std::pair<int, Vector3d>>,  std::vector<std::pair<int, Vector3d>>>> &manybodyNeighborIndex,
                                  const std::vector<std::pair<int, Vector3d>> &Ni,std::vector<std::pair<int, Vector3d>> &currentOriginalNeighbrs,  int order, bool saveBothWays,const int maxOrder);

    std::vector<std::pair<int, Vector3d>> getIntersection(const std::vector<std::pair<int, Vector3d>> &Ni, const std::vector<std::pair<int, Vector3d>> &Nj)
    {
        std::vector<std::pair<int, Vector3d>> N_intersection;
        N_intersection.reserve(Ni.size()/2);
        std::set_intersection(Ni.begin(), Ni.end(),
                              Nj.begin(), Nj.end(),
                              std::back_inserter(N_intersection),
                              NeighborPairCompare());
        return N_intersection;
    }

    std::vector<std::pair<int, Vector3d>> translateAllNi(std::vector<std::pair<int, Vector3d>> &Ni, const Vector3d &unitCellOffset) const ;



    // Lets put this on ice since it might be a bit dangerous to have this custom comparison
    // std::vector<std::pair<int, Vector3d>> getIntersectionOffset(const std::vector<std::pair<int, Vector3d>> &Ni, const std::vector<std::pair<int, Vector3d>> &Nj, const Vector3d &offset_j)
    // {
    //     std::vector<std::pair<int, Vector3d>> N_intersection;
    //     std::set_intersection(Ni.begin(), Ni.end(),
    //                           Nj.begin(), Nj.end(),
    //                           std::back_inserter(N_intersection),
    //                           NeighborPairCompareOffset(offset_j));
    //     return N_intersection;
    // }



  private:
    std::vector<double> _cutoffs;
};