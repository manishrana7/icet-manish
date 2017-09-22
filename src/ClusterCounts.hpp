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
#include "PeriodicTable.hpp"
#include "OrbitList.hpp"

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
  void countLatticeNeighbors(const Structure &,
                             const std::vector<std::pair<std::vector<LatticeNeighbor>, std::vector<LatticeNeighbor>>> &);
  void count(const Structure &,
             const std::vector<LatticeNeighbor> &);
  void count(const Structure &, const std::vector<std::vector<LatticeNeighbor>> &,
             const Cluster &);
  void countCluster(const Cluster &, const std::vector<int> &);
  void countOrbitlist(const Structure &, const OrbitList &);
  std::unordered_map<Cluster, std::map<std::vector<int>, int>> getClusterCounts() const
  {
    return _clusterCounts;
  }

  

  size_t size() const
  {
    return _clusterCounts.size();
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
      int total = 0;
      map_pair.first.print();
      std::cout << "==============" << std::endl;
      for (const auto &element_count_pair : map_pair.second)
      {
        for (auto el : element_count_pair.first)
        {
          std::cout << PeriodicTable::intStr[el] << " ";
        }
        std::cout << map_pair.second.at(element_count_pair.first) << std::endl;
        total +=map_pair.second.at(element_count_pair.first);
      }
      std::cout<<"Total: "<<total<<std::endl;
      std::cout << std::endl;
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