#include "Structure.hpp"
#include "Neighborlist.hpp"
#include "ManybodyNeighborlist.hpp"
#include "Cluster.hpp"
#include "PermutationMap.hpp"
#include "LatticeNeighbor.hpp"
#include "ClusterCounts.hpp"
#include <pybind11/pybind11.h>
#include "Orbit.hpp"
#include "OrbitList.hpp"
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <pybind11/operators.h>

PYBIND11_PLUGIN(_icetdev)
{
    py::module m("_icetdev", "pybind11 _icetdev plugin");

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
        .def("find_index_of_position_pybind", &Structure::findIndexOfPosition,
             py::arg("position"), py::arg("position_tolerance") = 1e-6)
        .def("findLatticeNeighborFromPosition", &Structure::findLatticeNeighborFromPosition,
             py::arg("position"), py::arg("position_tolerance") = 1e-6)
        .def("findLatticeNeighborsFromPositions", &Structure::findLatticeNeighborsFromPositions,
             py::arg("positions"), py::arg("position_tolerance") = 1e-6)
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
        .def("build", &ManybodyNeighborlist::build)
        ;

        

    py::class_<Cluster>(m, "Cluster")
        .def(py::init<std::vector<int> &, std::vector<double> &, const bool, const int>(), pybind11::arg("sites"),
             pybind11::arg("distances"), pybind11::arg("sortedCluster") = true, pybind11::arg("clusterTag") = 0)
        .def(py::init<const Structure &, const std::vector<LatticeNeighbor> &, const bool, const int>(), pybind11::arg("structure"),
             pybind11::arg("latticeNeighbors"), pybind11::arg("sortedCluster") = true, pybind11::arg("clusterTag") = 0)
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
        .def("print", &LatticeNeighbor::print )
        .def_readwrite("index", &LatticeNeighbor::index)
        .def_readwrite("unitcellOffset", &LatticeNeighbor::unitcellOffset)
        .def(py::self < py::self)        
        .def(py::self == py::self)
        .def(py::self + Eigen::Vector3d())
        ;

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

    py::class_<Orbit>(m, "Orbit")
        .def(py::init<const Cluster &>())
        .def("add_equivalent_sites",(void (Orbit::*)(const std::vector<LatticeNeighbor> &)) &Orbit::addEquivalentSites )
        .def("add_equivalent_sites",(void (Orbit::*)(const std::vector<std::vector<LatticeNeighbor>> &)) &Orbit::addEquivalentSites )
        .def("get_representative_cluster", &Orbit::getRepresentativeCluster)
        .def("get_equivalent_sites", &Orbit::getEquivalentSites)
        .def("size", &Orbit::size)
        .def("get_number_of_duplicates", &Orbit::getNumberOfDuplicates, py::arg("verbosity") = 0)
        .def(py::self < py::self)
        .def(py::self + Eigen::Vector3d())
        ;

    py::class_<OrbitList>(m, "OrbitList")
        .def(py::init<>())
        .def(py::init<const std::vector<Neighborlist> & , const Structure &>())
        .def(py::init<const Structure &, const std::vector<std::vector<LatticeNeighbor>> &, const std::vector<Neighborlist> &>())
        .def("add_orbit", &OrbitList::addOrbit)
        .def("get_number_of_NClusters", &OrbitList::getNumberOfNClusters)
        .def("get_orbit", &OrbitList::getOrbit)
        .def("clear", &OrbitList::clear)
        .def("sort", &OrbitList::sort)
        .def("get_orbitList", &OrbitList::getOrbitList)
        .def("size", &OrbitList::size)
        .def("print", &OrbitList::print, py::arg("verbosity") = 0)
        ;

    return m.ptr();
}
