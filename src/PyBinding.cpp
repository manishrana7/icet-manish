#include "Structure.hpp"
#include "Neighborlist.hpp"
#include "ManybodyNeighborlist.hpp"
#include "Cluster.hpp"
#include "PermutationMap.hpp"
#include "LatticeSite.hpp"
#include "ClusterCounts.hpp"
#include "LocalOrbitlistGenerator.hpp"
#include "ClusterSpace.hpp"
#include <pybind11/pybind11.h>
#include "Symmetry.hpp"
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
        .def(py::init<>())
        .def(py::init<const Eigen::Matrix<double, Dynamic, 3, Eigen::RowMajor> &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix3d &,
                      const std::vector<bool> &,
                      double>(),
             "Initialize a structure.",
             py::arg("positions"),
             py::arg("elements"),
             py::arg("cell"),
             py::arg("pbc"),
             py::arg("precision") = 1e-5)
        .def("get_positions", &Structure::getPositions)
        .def("set_positions", &Structure::setPositions)
        .def_property("positions",
                      &Structure::getPositions,
                      &Structure::setPositions,
                      "list of lists : atomic positions")
        .def("get_atomic_numbers", &Structure::getAtomicNumbers)
        .def("set_atomic_numbers", &Structure::setAtomicNumbers)
        .def_property("atomic_numbers",
                      &Structure::getAtomicNumbers,
                      &Structure::setAtomicNumbers,
                      "list of ints : atomic numbers of species on each site")
        .def("get_chemical_symbols", &Structure::getChemicalSymbols)
        .def("set_chemical_symbols", &Structure::setChemicalSymbols)
        .def_property("chemical_symbols",
                      &Structure::getChemicalSymbols,
                      &Structure::setChemicalSymbols,
                      "list of strings : chemical symbols of species on each site")
        .def("set_unique_sites", &Structure::setUniqueSites)
        .def("get_unique_sites", &Structure::getUniqueSites)
        .def_property("unique_sites",
                      &Structure::getUniqueSites,
                      &Structure::setUniqueSites,
                      "list of ints : unique sites")
        .def("get_unique_site", &Structure::getUniqueSite)
        .def("get_position", &Structure::getPosition)
        .def("get_distance", &Structure::getDistance)
        .def("get_distance2", &Structure::getDistance2)
        .def("find_site_by_position",
             &Structure::findSiteByPosition,
             py::arg("position"),
             "Returns the index of the site that matches the position.\n\n"
             "Parameters\n----------\n"
             "position : list/NumPy array\n"
             "    position in fractional coordinates")
        .def("find_lattice_neighbor_by_position",
             &Structure::findLatticeSiteByPosition,
             py::arg("position"))
        .def("find_lattice_neighbors_by_positions",
             &Structure::findLatticeSitesByPositions,
             py::arg("positions"))
        .def("get_pbc", &Structure::getPBC)
        .def("set_pbc", &Structure::setPBC)
        .def_property("pbc",
                      &Structure::getPBC,
                      &Structure::setPBC,
                      "3-dimensional vector : periodic boundary conditions")
        .def("get_cell", &Structure::getCell)
        .def("set_cell", &Structure::setCell)
        .def_property("cell",
                      &Structure::getCell,
                      &Structure::setCell,
                      "3x3 array : cell metric")
        .def("__len__", &Structure::size);

    py::class_<Neighborlist>(m, "Neighborlist")
        .def(py::init<const double>())
        .def("build", &Neighborlist::build)
        .def("is_neighbor", &Neighborlist::isNeighbor)
        .def("get_neighbors", &Neighborlist::getNeighbors)
        .def("__len__", &Neighborlist::size)

        ;

    py::class_<ManybodyNeighborlist>(m, "ManybodyNeighborlist")
        .def(py::init<>())
        .def("calc_intersection", &ManybodyNeighborlist::getIntersection)
        .def("build", &ManybodyNeighborlist::build);

    py::class_<Cluster>(m, "Cluster")
        // .def(py::init<std::vector<int> &, std::vector<double> &, const bool, const int>(), pybind11::arg("sites"),
        //      pybind11::arg("distances"), pybind11::arg("sortedCluster") = true, pybind11::arg("clusterTag") = 0)
        .def(py::init<const Structure &, const std::vector<LatticeSite> &, const bool, const int>(), pybind11::arg("structure"),
             pybind11::arg("latticeNeighbors"), pybind11::arg("sortedCluster") = true, pybind11::arg("clusterTag") = 0)
        .def("count", &Cluster::count)
        .def("get_count", &Cluster::getCount)
        .def("get_sites", &Cluster::getSites)
        .def("get_distances", &Cluster::getDistances)
        .def("print", &Cluster::print)
        .def("is_sorted", &Cluster::isSorted)
        .def("get_clustertag", &Cluster::getClusterTag)
        .def("get_geometrical_size", &Cluster::getGeometricalSize)
        .def("get_number_of_bodies",&Cluster::getNumberOfBodies)
        .def("__hash__", [](const Cluster &cluster) { return std::hash<Cluster>{}(cluster); })
        .def(py::self < py::self)
        .def(py::self == py::self)
        // .def(hash(py::self))
        ;

    py::class_<PermutationMap>(m, "PermutationMap")
        .def(py::init<const std::vector<Vector3d> &,
                      const std::vector<Matrix3d> &>())
        .def("build", &PermutationMap::build)
        .def("get_permutated_positions", &PermutationMap::getPermutatedPositions)
        .def("get_indiced_positions", &PermutationMap::getIndicedPermutatedPositions)

        ;

    py::class_<LatticeSite>(m, "LatticeSite")
        .def(py::init<const int, const Vector3d &>())
        .def("print", &LatticeSite::print)
        .def_property("index", &LatticeSite::index, &LatticeSite::setIndex)
        .def_property("unitcellOffset", &LatticeSite::unitcellOffset, &LatticeSite::setUnitcellOffset)
        .def(py::self < py::self)
        .def(py::self == py::self)
        .def(py::self + Eigen::Vector3d())
        .def("__hash__", [](const LatticeSite &latticeNeighbor) { return std::hash<LatticeSite>{}(latticeNeighbor); })

        ;

    py::class_<ClusterCounts>(m, "ClusterCounts")
        .def(py::init<>())
        .def("count_lattice_neighbors", &ClusterCounts::countLatticeSites)
        .def("count", (void (ClusterCounts::*)(const Structure &, const std::vector<LatticeSite> &)) & ClusterCounts::count)
        .def("count", (void (ClusterCounts::*)(const Structure &, const std::vector<std::vector<LatticeSite>> &, const Cluster &)) & ClusterCounts::count)
        .def("count_orbitlist",&ClusterCounts::countOrbitlist)
        .def("__len__", &ClusterCounts::size)
        .def("reset", &ClusterCounts::reset)
        .def("get_cluster_counts", [](const ClusterCounts &clusterCounts) {
            //&ClusterCounts::getClusterCounts
            py::dict clusterCountDict;
            for (const auto &mapPair : clusterCounts.getClusterCounts())
            {
                py::dict d;
                for (const auto &vecInt_intPair : mapPair.second)
                {
                    d[py::tuple(py::cast(vecInt_intPair.first))] = vecInt_intPair.second;
                }
                clusterCountDict[py::cast(mapPair.first)] = d;
            }
            return clusterCountDict;
        })
        .def("print", &ClusterCounts::print);

    py::class_<Orbit>(m, "Orbit")
        .def(py::init<const Cluster &>())
        .def("add_equivalent_sites", (void (Orbit::*)(const std::vector<LatticeSite> &, bool)) & Orbit::addEquivalentSites, py::arg("lattice_neighbors"), py::arg("sort")=false)
        .def("add_equivalent_sites", (void (Orbit::*)(const std::vector<std::vector<LatticeSite>> &, bool)) & Orbit::addEquivalentSites, py::arg("lattice_neighbors"), py::arg("sort")=false)
        .def("get_representative_cluster", &Orbit::getRepresentativeCluster)
        .def("get_equivalent_sites", &Orbit::getEquivalentSites)
        .def("get_representative_sites", &Orbit::getRepresentativeSites)
        .def("get_equivalent_sites_permutations", &Orbit::getEquivalentSitesPermutations)
        .def("get_sites_with_permutation", &Orbit::getSitesWithPermutation)
        .def("__len__", &Orbit::size)
        .def("get_number_of_duplicates", &Orbit::getNumberOfDuplicates, py::arg("verbosity") = 0)
        .def("get_MC_vectors",&Orbit::getMCVectors)
        .def(py::self < py::self)
        .def(py::self == py::self)
        .def(py::self + Eigen::Vector3d());

    py::class_<OrbitList>(m, "OrbitList")
        .def(py::init<>())
        .def(py::init<const std::vector<Neighborlist> &, const Structure &>())
        .def(py::init<const Structure &, const std::vector<std::vector<LatticeSite>> &, const std::vector<Neighborlist> &>())
        .def("add_orbit", &OrbitList::addOrbit)
        .def("get_number_of_NClusters", &OrbitList::getNumberOfNClusters)
        .def("get_orbit", &OrbitList::getOrbit)
        .def("clear", &OrbitList::clear)
        .def("sort", &OrbitList::sort)
        .def("get_orbitList", &OrbitList::getOrbitList)
        .def("get_primitive_structure",&OrbitList::getPrimitiveStructure)
        .def("__len__", &OrbitList::size)
        .def("print", &OrbitList::print, py::arg("verbosity") = 0)
        .def("get_supercell_orbitlist", &OrbitList::getSupercellOrbitlist);

    py::class_<LocalOrbitlistGenerator>(m, "LocalOrbitlistGenerator")
        .def(py::init<const OrbitList &, const Structure &>())
        .def("generate_local_orbitlist", (OrbitList(LocalOrbitlistGenerator::*)(const unsigned int)) & LocalOrbitlistGenerator::generateLocalOrbitlist)
        .def("generate_local_orbitlist", (OrbitList(LocalOrbitlistGenerator::*)(const Vector3d &)) & LocalOrbitlistGenerator::generateLocalOrbitlist)
        .def("clear", &LocalOrbitlistGenerator::clear)
        .def("get_unique_offsets_count", &LocalOrbitlistGenerator::getUniqueOffsetsCount)
        .def("get_prim_to_supercell_map", &LocalOrbitlistGenerator::getPrimToSupercellMap)
        .def("get_unique_primcell_offsets", &LocalOrbitlistGenerator::getUniquePrimcellOffsets);

        py::class_<ClusterSpace>(m, "ClusterSpace",py::dynamic_attr())
        .def(py::init<std::vector<int>, std::vector<std::string>, const OrbitList &>())
        .def("_get_clustervector",&ClusterSpace::generateClustervector)
        .def("get_orbit", &ClusterSpace::getOrbit)
        .def("get_cluster_product", &ClusterSpace::getClusterProduct)
        .def("get_clusterspace_info", &ClusterSpace::getClusterSpaceInfo)
        .def("get_clusterspace_size", &ClusterSpace::getClusterSpaceSize)
        .def("get_elements", &ClusterSpace::getAtomicNumbers)
        .def("get_cutoffs",&ClusterSpace::getCutoffs)
        .def("get_primitive_structure",&ClusterSpace::getPrimitiveStructure)
        .def("get_native_clusters",&ClusterSpace::getNativeClusters)
        ;


    return m.ptr();
}
