#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include "Vector3dCompare.hpp"
#include "Neighborlist.hpp"
#include <vector>
#include "LatticeNeighbor.hpp"
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

    std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> build(const std::vector<Neighborlist> &, int index, bool);

    void combineToHigherOrder(const Neighborlist &nl,
                              std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &manybodyNeighborIndex,
                              const std::vector<LatticeNeighbor> &Ni, std::vector<LatticeNeighbor> &currentOriginalNeighbrs,  bool saveBothWays, const int);

    std::vector<LatticeNeighbor> getIntersection(const std::vector<LatticeNeighbor> &Ni, const std::vector<LatticeNeighbor> &Nj)
    {
        std::vector<LatticeNeighbor> N_intersection;
        N_intersection.reserve(Ni.size());
        std::set_intersection(Ni.begin(), Ni.end(),
                              Nj.begin(), Nj.end(),
                              std::back_inserter(N_intersection));
        return N_intersection;
    }

    void translateAllNi(std::vector<LatticeNeighbor> &Ni, const Vector3d &unitCellOffset) const;

  private:
    std::vector<double> _cutoffs;
    std::vector<LatticeNeighbor> getFilteredNj(const std::vector<LatticeNeighbor> &, const LatticeNeighbor &) const;
};