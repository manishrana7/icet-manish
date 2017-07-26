#include "Structure.hpp"
#include "Neighborlist.hpp"
#include "ManybodyNeighborlist.hpp"
#include "Cluster.hpp"
#include "PermutationMap.hpp"
#include "LatticeNeighbor.hpp"
#include "ClusterCounts.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <pybind11/operators.h>

PYBIND11_PLUGIN(example)
{
  py::module m("example", "pybind11 example plugin");

  py::class_<Structure>(m, "Structure")
      .def(py::init<const Eigen::Matrix<double, Dynamic, 3, Eigen::RowMajor> &,
                    const std::vector<std::string> &,
                    const Eigen::Matrix3d &,
                    const std::vector<bool> &>())
      .def("set_positions", &Structure::setPositions)
      .def("set_elements", &Structure::setElements)
      .def("get_elements", &Structure::getElements)
      .def("set_unique_sites", &Structure::setUniqueSites)
      .def("get_unique_sites", &Structure::getUniqueSites)
      .def("get_unique_site", &Structure::getSite)
      .def("get_positions", &Structure::getPositions)
      .def("get_position", &Structure::getPosition)
      .def("get_distance", &Structure::getDistance)
      .def("get_distance2", &Structure::getDistance2)
      .def("has_pbc", &Structure::has_pbc)
      .def("get_pbc", &Structure::get_pbc)
      .def("set_pbc", &Structure::set_pbc)
      .def("get_cell", &Structure::get_cell)
      .def("set_cell", &Structure::set_cell)
      .def("size", &Structure::size);

  py::class_<Neighborlist>(m, "Neighborlist")
      .def(py::init<const double>())
      .def("build", &Neighborlist::build)
      .def("is_neighbor", &Neighborlist::isNeighbor)
      .def("get_neighbors", &Neighborlist::getNeighbors)
      .def("size", &Neighborlist::size)

      ;

  py::class_<ManybodyNeighborlist>(m, "ManybodyNeighborlist")
      .def(py::init<>())
      .def("calc_intersection", &ManybodyNeighborlist::getIntersection)
      .def("build", &ManybodyNeighborlist::build);

  py::class_<Cluster>(m, "Cluster")
      .def(py::init<std::vector<int> &, std::vector<double> &, const bool, const int>(), pybind11::arg("sites"),
           pybind11::arg("distances"), pybind11::arg("sortedCluster") = true, pybind11::arg("clusterTag") = 0)
      .def("count", &Cluster::count)
      .def("get_count", &Cluster::getCount)
      .def("get_sites", &Cluster::getSites)
      .def("get_distances", &Cluster::getDistances)
      .def("print", &Cluster::print)
      .def("is_sorted", &Cluster::isSorted)
      .def("get_clusterTag", &Cluster::getClusterTag)
      .def(py::self < py::self);

  py::class_<PermutationMap>(m, "PermutationMap")
      .def(py::init<const std::vector<Vector3d> &,
                    const std::vector<Matrix3d> &>())
      .def("build", &PermutationMap::build)
      .def("get_permutated_positions", &PermutationMap::getPermutatedPositions)
      .def("get_indiced_positions", &PermutationMap::getIndicedPermutatedPositions)

      ;

  py::class_<LatticeNeighbor>(m, "LatticeNeighbor")
      .def(py::init<const int, const Vector3d &>())
      .def_readwrite("index", &LatticeNeighbor::index)
      .def_readwrite("unitcellOffset", &LatticeNeighbor::unitcellOffset)
      .def(py::self < py::self);

  py::class_<ClusterCounts>(m, "ClusterCounts")
      .def(py::init<>())
      .def("count_lattice_neighbors", &ClusterCounts::countLatticeNeighbors)
      .def("count_singlets", &ClusterCounts::countSinglets)
      .def("count_pairs", &ClusterCounts::countPairs)
      .def("count", &ClusterCounts::count)
      .def("size", &ClusterCounts::size)
      .def("reset", &ClusterCounts::reset)
      .def("get_cluster_counts", &ClusterCounts::getClusterCounts)
      .def("print", &ClusterCounts::print);

  return m.ptr();
}
