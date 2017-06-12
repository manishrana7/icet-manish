#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <unordered_set>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "ManybodyNeighborlist.hpp"
#include "Structure.hpp"
#include "Cluster.hpp"
#include "LatticeNeighbor.hpp"

using namespace Eigen;

namespace py = pybind11;

class ClusterCounts
{
public:
  ClusterCounts()
  {
    symprec = 1e-5;
    //empty constructor
  }

  //void count(const Structure & , XXX indices );
  void count_using_mbnl(const Structure &, ManybodyNeighborlist &, const int);
  void countLatticeNeighbors(const Structure &,
                             const std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &);
  void count(const Structure &,
             const std::vector<LatticeNeighbor> &);
  std::unordered_map<Cluster, std::map<std::vector<int>, int>> getClusterCounts() const
  {
    return _clusterCounts;
  }
  void reset()
  {
    _clusterCounts.clear();
  }

  /**
 Helpful function that prints the cluster counts
*/
  void print()
  {
    for (const auto &map_pair : _clusterCounts)
    {
      map_pair.first.print();
      std::cout << "==============" << std::endl;
      for (const auto &element_count_pair : map_pair.second)
      {
        for (auto el : element_count_pair.first)
        {
          std::cout << el << " ";
        }
        std::cout << map_pair.second.at(element_count_pair.first) << std::endl;
      }
    }
  }

private:
  std::unordered_map<Cluster, std::map<std::vector<int>, int>> _clusterCounts;

  double roundDouble(const double &double_value)
  {
    return round(double_value * 1.0 / symprec) / (1.0 / symprec);
  }
  double symprec;
};