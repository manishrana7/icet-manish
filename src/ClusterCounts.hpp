#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "ManybodyNeighborlist.hpp"
#include "Structure.hpp"
#include "Cluster.hpp"
#include "unordered_map"
using namespace Eigen;

namespace py = pybind11;

class ClusterCounts
{
  public:
    ClusterCounts();

    //void count(const Structure & , XXX indices );
    void count_using_mbnl(const Structure &, ManybodyNeighborlist &, const int);

  private:
    void countLatticeNeighbors(const Structure &,
                               const std::vector<std::pair<std::vector<std::pair<int, Vector3d>>, std::vector<std::pair<int, Vector3d>>>> &);
    void count(const Structure &,
               const std::vector<std::pair<int, Vector3d>> &);
};