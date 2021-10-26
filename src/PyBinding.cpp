#include <iostream>
#include <sstream>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>

#include "Cluster.hpp"
#include "ClusterExpansionCalculator.hpp"
#include "ClusterSpace.hpp"
#include "LatticeSite.hpp"
#include "LocalOrbitListGenerator.hpp"
#include "ManyBodyNeighborList.hpp"
#include "Orbit.hpp"
#include "OrbitList.hpp"
#include "PeriodicTable.hpp"
#include "Structure.hpp"
#include "Symmetry.hpp"

PYBIND11_MODULE(_icet, m)
{

    m.doc() = R"pbdoc(
        Python-C++ interface
        ====================

        This is the Python interface generated via pybind11 from the C++
        core classes and methods.

        .. toctree::
           :maxdepth: 2

        .. currentmodule:: _icet

        Cluster
        -------
        .. autoclass:: Cluster
           :members:

        ClusterSpace
        ------------
        .. autoclass:: ClusterSpace
           :members:

        LatticeSite
        -----------
        .. autoclass:: LatticeSite
           :members:

        LocalOrbitListGenerator
        -----------------------
        .. autoclass:: LocalOrbitListGenerator
           :members:

        ManyBodyNeighborList
        --------------------
        .. autoclass:: ManyBodyNeighborList
           :members:

        Orbit
        -----
        .. autoclass:: Orbit
           :members:

        _OrbitList
        ----------
        .. autoclass:: _OrbitList
           :members:

        Structure
        ---------
        .. autoclass:: Structure
           :members:
    )pbdoc";

    // Disable the automatically generated signatures that prepend the
    // docstrings by default.
    py::options options;
    options.disable_function_signatures();

    py::class_<Structure>(m, "_Structure")
        .def(py::init<>())
        .def(py::init<const Eigen::Matrix<double, Dynamic, 3, Eigen::RowMajor> &,
                      const py::array_t<int> &,
                      const Eigen::Matrix3d &,
                      const std::vector<bool> &>(),
             "Initializes an icet Structure instance.",
             py::arg("positions"),
             py::arg("atomic_numbers"),
             py::arg("cell"),
             py::arg("pbc"))
        .def_property(
            "pbc",
            &Structure::getPBC,
            &Structure::setPBC,
            "list(int) : periodic boundary conditions")
        .def_property(
            "cell",
            &Structure::getCell,
            &Structure::setCell,
            "list(list(float)) : cell metric")
        .def_property(
            "positions",
            &Structure::getPositions,
            &Structure::setPositions,
            "list(list(float)) : atomic positions in Cartesian coordinates")
        .def_property("atomic_numbers",
                      &Structure::getAtomicNumbers,
                      &Structure::setAtomicNumbers,
                      "list(int) : atomic numbers of species on each site")
        .def("set_number_of_allowed_species",
             (void (Structure::*)(const std::vector<int> &)) & Structure::setNumberOfAllowedSpecies,
             py::arg("numbersOfAllowedSpecies"),
             R"pbdoc(
             Sets the number of allowed species on each site.

             This method allows one to specify for each site in the structure
             the number of species allowed on that site.

             Parameters
             ----------
             numbersOfAllowedSpecies : list(int)
             )pbdoc")
        .def("get_position",
             &Structure::getPosition,
             py::arg("site"),
             R"pbdoc(
             Returns the position of a specified site

             Parameters
             ----------
             site : LatticeSite object
                site of interest

             Returns
             -------
             vector
                 position in Cartesian coordinates
             )pbdoc")
        .def("get_distance",
             &Structure::getDistance,
             py::arg("index1"),
             py::arg("index2"),
             py::arg("offset1") = Vector3d(0, 0, 0),
             py::arg("offset2") = Vector3d(0, 0, 0),
             R"pbdoc(
             Returns the distance between two sites

             Parameters
             ----------
             index1 : int
                 index of the first site
             index2 : int
                 index of the second site
             offset1 : vector
                 offset to be applied to the first site
             offset2 : vector
                 offset to be applied to the second site

             Returns
             -------
             float
                 distance in length units
             )pbdoc")
        .def("find_lattice_site_by_position",
             &Structure::findLatticeSiteByPosition,
             R"pbdoc(
             Returns the lattice site that matches the position.

             Parameters
             ----------
             position : list or ndarray
                 position in Cartesian coordinates
             fractional_position_tolerance : float
                 tolerance for positions in fractional coordinates

             Returns
             -------
             _icet.LatticeSite
                 lattice site
             )pbdoc",
             py::arg("position"),
             py::arg("fractional_position_tolerance"))
        .def("__len__", &Structure::size);

    // @todo document ManyBodyNeighborList in pybindings
    py::class_<ManyBodyNeighborList>(m, "ManyBodyNeighborList",
                                     R"pbdoc(
        This class handles a many-body neighbor list.
        )pbdoc")
        .def(py::init<>())
        .def("calculate_intersection", &ManyBodyNeighborList::getIntersection)
        .def("build", &ManyBodyNeighborList::build);

    py::class_<Cluster>(m, "Cluster",
                        R"pbdoc(
        This class handles a many-body neighbor list.

        Parameters
        ----------
        structure : icet Structure instance
            atomic configuration
        lattice_sites : list(int)
            list of lattice sites that form the cluster
        )pbdoc")
        .def(py::init<const Structure &,
                      const std::vector<LatticeSite> &>(),
             "Initializes a cluster instance.",
             py::arg("structure"),
             py::arg("lattice_sites"))
        .def_property_readonly(
            "distances",
            &Cluster::distances,
            "list(float) : list of distances between sites")
        .def_property_readonly(
            "radius",
            &Cluster::radius,
            "float : the radius of the cluster")
        .def_property_readonly(
            "order",
            &Cluster::order,
            "int : order of the cluster (= number of sites)")
        .def("__len__",
             &Cluster::order)
        .def("__str__",
             [](const Cluster &cluster)
             {
                 std::ostringstream msg;
                 msg << "radius: " << cluster.radius();
                 msg << " vertices:";
                 for (const auto dist : cluster.distances())
                 {
                     msg << " " << std::to_string(dist);
                 }
                 return msg.str();
             });
    ;

    py::class_<LatticeSite>(m, "LatticeSite",
                            R"pbdoc(
        This class handles a lattice site.

        Parameters
        ----------

        )pbdoc")
        .def(py::init<const int,
                      const Vector3d &>(),
             "Initializes a LatticeSite object.",
             py::arg("site_index"),
             py::arg("unitcell_offset"))
        .def_property(
            "index",
            &LatticeSite::index,
            &LatticeSite::setIndex,
            "int : site index")
        .def_property(
            "unitcell_offset",
            &LatticeSite::unitcellOffset,
            &LatticeSite::setUnitcellOffset,
            "list(int) : unit cell offset (in units of the cell vectors)")
        .def(py::self < py::self)
        .def(py::self == py::self)
        .def(py::self + Eigen::Vector3d())
        .def("__hash__", [](const LatticeSite &latticeNeighbor)
             { return std::hash<LatticeSite>{}(latticeNeighbor); });

    // @todo convert getters to properties
    // @todo document Orbit in pybindings
    py::class_<Orbit>(m, "Orbit")
        .def(py::init<const Cluster &>())
        .def_property_readonly(
            "representative_cluster",
            &Orbit::getRepresentativeCluster,
            "cluster to which all other symmetry equivalent clusters can be related")
        .def_property_readonly(
            "sites_of_representative_cluster", &Orbit::getSitesOfRepresentativeCluster,
            "list of sites that comprise the representative cluster")
        .def_property_readonly(
            "order",
            [](const Orbit &orbit)
            { return orbit.order(); },
            "number of sites in the representative cluster")
        .def_property_readonly(
            "radius",
            [](const Orbit &orbit)
            { return orbit.radius(); },
            "radius of the representative cluster")
        .def_property(
            "permutations_to_representative",
            &Orbit::getPermutationsOfEquivalentClusters,
            &Orbit::setPermutationsOfEquivalentClusters,
            R"pbdoc(
             list of permutations;
             permutations_to_representative[i] takes self.equivalent_clusters[i] to
             the same sorting as self.representative_cluster.

             This can be used if you for example want to count elements and are
             interested in difference between ABB, BAB, BBA and so on. If you count
             the lattice sites that are permuted according to these permutations
             then you will get the correct counts.
             )pbdoc")
        .def_property(
            "allowed_permutations",
            [](const Orbit &orbit)
            {
                std::set<std::vector<int>> allowedPermutations = orbit.getAllowedClusterPermutations();
                std::vector<std::vector<int>> retPermutations(allowedPermutations.begin(), allowedPermutations.end());
                return retPermutations;
            },
            [](Orbit &orbit, const std::vector<std::vector<int>> &newPermutations)
            {
                std::set<std::vector<int>> allowedPermutations;
                for (const auto &perm : newPermutations)
                {
                    allowedPermutations.insert(perm);
                }
                orbit.setAllowedClusterPermutations(allowedPermutations);
            },
            R"pbdoc(
             Gets the list of equivalent permutations for this orbit. If this
             orbit is a triplet and the permutation [0,2,1] exists this means
             that The lattice sites [s1, s2, s3] are equivalent to [s1, s3,
             s2] This will have the effect that for a ternary CE the cluster
             functions (0,1,0) will not be considered since it is equivalent
             to (0,0,1).
             )pbdoc")
        .def_property(
            "equivalent_clusters",
            &Orbit::getEquivalentClusters,
            &Orbit::setEquivalentClusters,
            "list of symmetry equivalent clusters")
        .def_property_readonly(
            "permuted_equivalent_clusters",
            &Orbit::getPermutedEquivalentClusters,
            "equivalent clusters permuted to match the sorting of the representative cluster")
        .def("get_permuted_cluster_by_index",
             &Orbit::getPermutedClusterByIndex,
             R"pbdoc(
             Returns the equivalent cluster at position `index` using
             the permutation of the representative cluster.

             Parameters
             ----------
             index : int
                index of site to return
             )pbdoc",
             py::arg("index"))
        .def("get_mc_vectors", &Orbit::getMultiComponentVectors,
             R"pbdoc(
             Return the multi-component vectors for this orbit given the allowed components.
             The multi-component vectors are returned as a list of tuples.

             Parameters
             ----------
             allowed_components : list(int)
                The allowed components for the lattice sites,
                allowed_components[i] correspond to the number
                of allowed compoments at lattice site
                orbit.representative_cluster[i].)pbdoc")
        .def(
            "count_clusters",
            [](const Orbit &orbit,
               const Structure &structure,
               const int siteIndexForDoubleCountingCorrection,
               const int permuteClusters)
            {
                py::dict clusterCountDict;
                for (const auto &mapPair : orbit.countClusters(structure,
                                                               siteIndexForDoubleCountingCorrection,
                                                               permuteClusters))
                {
                    py::list element_symbols;
                    for (auto el : mapPair.first)
                    {
                        auto getElementSymbols = PeriodicTable::intStr[el];
                        element_symbols.append(getElementSymbols);
                    }
                    double countDouble = mapPair.second;
                    if (std::abs(std::round(countDouble) - countDouble) > 1e-6)
                    {
                        std::runtime_error("Cluster count is a non-integer.");
                    }
                    int count = (int)std::round(countDouble);
                    clusterCountDict[py::tuple(element_symbols)] = count;
                }
                return clusterCountDict;
            },
            R"pbdoc(
             Count clusters in this orbit for a structure.

             Parameters
             ----------
             structure : Structure
                Structure to count clusters for
             site_index_for_double_counting_correction : int
                Avoid double counting clusters containing this index
                (default -1, no such correction)
             permute_clusters : bool
                Permute clusters before counting (default: false)
             )pbdoc",
            py::arg("structure"),
            py::arg("site_index_for_double_counting_correction") = -1,
            py::arg("permute_sites") = false)
        .def("sort", &Orbit::sort,
             "Sorts the list of equivalent sites.")
        .def("get_all_possible_mc_vector_permutations",
             &Orbit::getAllPossibleMultiComponentVectorPermutations,
             R"pbdoc(
             Similar to get all permutations but needs to be filtered through the number of allowed elements.

             Parameters
             ----------
             allowed_components : list(int)
                 The allowed components for the lattice sites,
                 `allowed_components[i]` correspond to the lattice site
                 `self.representative_cluster[i]`.

             returns all_mc_vectors : list(list(int)
             )pbdoc")
        .def("__len__", &Orbit::size)
        .def("__str__",
             [](const Orbit &orbit)
             {
                 std::ostringstream msg;
                 msg << "order: " << orbit.order() << std::endl;
                 msg << "multiplicity: " << orbit.size() << std::endl;
                 msg << "radius: " << orbit.radius() << std::endl;
                 msg << "representative_cluster:" << std::endl;
                 for (const auto site : orbit.getSitesOfRepresentativeCluster())
                 {
                     msg << "    site: " << site << std::endl;
                 }
                 msg << "equivalent_clusters:" << std::endl;
                 int k = -1;
                 for (const auto sites : orbit.getEquivalentClusters())
                 {
                     k += 1;
                     msg << "  cluster: " << k << std::endl;
                     for (const auto site : sites)
                     {
                         msg << "    site: " << site << std::endl;
                     }
                 }
                 return msg.str();
             })
        .def(py::self < py::self)
        .def(py::self + Eigen::Vector3d())
        .def(py::self += py::self);

    py::class_<OrbitList>(m, "_OrbitList",
                          R"pbdoc(
        This class manages an orbit list. The orbit list is constructed for the given
        structure using the matrix of equivalent sites and a list of neighbor lists.

        Parameters
        ----------
        structure : _icet.Structure
            (supercell) structure for which to generate orbit list
        matrix_of_equivalent_sites : list(list(_icet.LatticeSite))
            matrix of symmetry equivalent sites
        neighbor_lists : list(list(list(_icet.LatticeSite)))
            neighbor lists for each (cluster) order
        position_tolerance
            tolerance applied when comparing positions in Cartesian coordinates
        )pbdoc")
        .def(py::init<>())
        .def(py::init<const Structure &,
                      const std::vector<std::vector<LatticeSite>> &,
                      const std::vector<std::vector<std::vector<LatticeSite>>> &,
                      const double>(),
             "Constructs an OrbitList object from a matrix of equivalent sites.",
             py::arg("structure"),
             py::arg("matrix_of_equivalent_sites"),
             py::arg("neighbor_lists"),
             py::arg("position_tolerance"))
        .def_property_readonly(
            "orbits",
            &OrbitList::getOrbits,
            "list(_icet.Orbit) : list of orbits")
        .def("get_orbit_list", &OrbitList::getOrbits,
             "Returns the list of orbits")
        .def("add_orbit",
             &OrbitList::addOrbit,
             "Adds an orbit.")
        .def("get_orbit",
             &OrbitList::getOrbit,
             "Returns the orbit at position i in the orbit list.")
        .def("_remove_inactive_orbits",
             &OrbitList::removeInactiveOrbits)
        .def("clear",
             &OrbitList::clear,
             "Clears the list of orbits.")
        .def("sort", &OrbitList::sort,
             R"pbdoc(
             Sorts the orbits by order and radius.

             Parameters
             ----------
             position_tolerance : float
                 tolerance applied when comparing positions in Cartesian coordinates
             )pbdoc",
             py::arg("position_tolerance"))
        .def("remove_orbit",
             &OrbitList::removeOrbit,
             R"pbdoc(
             Removes the orbit with the input index.

             Parameters
             ---------
             index : int
                 index of the orbit to be removed
             )pbdoc")
        .def("_is_row_taken",
             &OrbitList::isRowsTaken,
             R"pbdoc(
             Returns true if rows exist in taken_rows.

             Parameters
             ----------
             taken_rows : set(tuple(int))
                 unique collection of row index
             rows : list(int)
                 row indices
             )pbdoc",
             py::arg("taken_rows"),
             py::arg("rows"))
        .def("_get_sites_translated_to_unitcell",
             &OrbitList::getSitesTranslatedToUnitcell,
             R"pbdoc(
             Returns a set of sites where at least one site is translated inside the unit cell.

             Parameters
             ----------
             lattice_neighbors : list(_icet.LatticeSite)
                set of lattice sites that might be representative for a cluster
             sort : bool
                If true sort translasted sites.
             )pbdoc",
             py::arg("lattice_neighbors"),
             py::arg("sort"))
        .def("_get_all_columns_from_sites",
             &OrbitList::getAllColumnsFromCluster,
             R"pbdoc(
             Finds the sites in column1, extract and return all columns along with their unit cell
             translated indistinguishable sites.

             Parameters
             ----------
             sites : list(_icet.LatticeSite)
                 sites that correspond to the columns that will be returned
             )pbdoc",
             py::arg("sites"))
        .def("get_primitive_structure",
             &OrbitList::getPrimitiveStructure,
             "Returns the primitive atomic structure used to construct the OrbitList instance.")
        .def("__len__",
             &OrbitList::size,
             "Returns the total number of orbits counted in the OrbitList instance.")
        .def_property_readonly("matrix_of_equivalent_positions",
                               &OrbitList::getMatrixOfEquivalentSites,
                               "list(list(_icet.LatticeSite)) : matrix_of_equivalent_positions");

    py::class_<LocalOrbitListGenerator>(m, "LocalOrbitListGenerator",
                                        R"pbdoc(
        This class handles the generation of local orbit lists, which are used in
        the computation of cluster vectors of supercells of the primitive structure.
        Upon initialization a LocalOrbitListGenerator object is constructed from an
        orbit list and a supercell structure.

        Parameters
        ----------
        orbit_list : _icet.OrbitList
            an orbit list set up from a primitive structure
        structure : _icet.Structure
            supercell build up from the same primitive structure used to set the input orbit list
        fractional_position_tolerance : float
            tolerance for positions in fractional coordinates
        )pbdoc")
        .def(py::init<const OrbitList &,
                      const Structure &,
                      const double>(),
             "Constructs a LocalOrbitListGenerator object from an orbit list and a structure.",
             py::arg("orbit_list"),
             py::arg("structure"),
             py::arg("fractional_position_tolerance"))
        .def("generate_local_orbit_list",
             (OrbitList(LocalOrbitListGenerator::*)(const size_t)) & LocalOrbitListGenerator::getLocalOrbitList,
             R"pbdoc(
             Generates and returns the local orbit list from an input index corresponding a specific offset of
             the primitive structure.

             Parameters
             ----------
             index : int
                 index of the unique offsets list
             )pbdoc",
             py::arg("index"))
        .def("generate_full_orbit_list",
             &LocalOrbitListGenerator::getFullOrbitList,
             R"pbdoc(
             Generates and returns a local orbit list, which orbits included the equivalent sites
             of all local orbit list in the supercell.
             )pbdoc")
        .def("get_number_of_unique_offsets",
             &LocalOrbitListGenerator::getNumberOfUniqueOffsets,
             "Returns the number of unique offsets")
        .def("_get_primitive_to_supercell_map",
             &LocalOrbitListGenerator::getMapFromPrimitiveToSupercell,
             "Returns the primitive to supercell mapping")
        .def("_get_unique_primcell_offsets",
             &LocalOrbitListGenerator::getUniquePrimitiveCellOffsets,
             "Returns a list with offsets of primitive structure that span to position of atoms in the supercell.");

    /// @todo Check which of the following members must actually be exposed.
    /// @todo Turn getters into properties if possible. (Some might require massaging in cluster_space.py.)
    py::class_<ClusterSpace>(m, "ClusterSpace", py::dynamic_attr())
        .def(py::init<std::vector<std::vector<std::string>> &,
                      const OrbitList,
                      const double,
                      const double>(),
             "Initializes an icet ClusterSpace instance.",
             py::arg("chemical_symbols"),
             py::arg("orbit_list"),
             py::arg("position_tolerance"),
             py::arg("fractional_position_tolerance"))
        .def(
            "get_cluster_vector",
            [](const ClusterSpace &clusterSpace,
               const Structure &structure,
               const double fractionalPositionTolerance)
            {
                auto cv = clusterSpace.getClusterVector(structure, fractionalPositionTolerance);
                return py::array(cv.size(), cv.data());
            },
            R"pbdoc(
             Returns the cluster vector corresponding to the input structure.
             The first element in the cluster vector will always be one (1) corresponding to
             the zerolet. The remaining elements of the cluster vector represent averages
             over orbits (symmetry equivalent clusters) of increasing order and size.

             Parameters
             ----------
             structure : _icet.Structure
                 atomic configuration
             fractional_position_tolerance : float
                 tolerance applied when comparing positions in fractional coordinates
             )pbdoc",
            py::arg("structure"),
            py::arg("fractional_position_tolerance"))
        .def(
            "_merge_orbit",
            &ClusterSpace::mergeOrbits,
            R"pbdoc(
             Merges two orbits. This implies that the equivalent clusters
             from the second orbit are added to to the list of equivalent
             clusters of the first orbit, after which the second orbit is
             removed.

             Parameters
             ----------
             index1 : int
                 index of the first orbit in the orbit list of the cluster space
             index2 : int
                 index of the second orbit in the orbit list of the cluster space
             )pbdoc",
            py::arg("index1"),
            py::arg("index2"))

        .def("_get_orbit_list", &ClusterSpace::getPrimitiveOrbitList)
        .def("get_orbit", &ClusterSpace::getOrbit)
        .def_property_readonly("species_maps", &ClusterSpace::getSpeciesMaps)
        .def("get_multi_component_vectors_by_orbit", &ClusterSpace::getMultiComponentVectorsByOrbit)
        .def("get_chemical_symbols",
             &ClusterSpace::getChemicalSymbols,
             "Returns list of species associated with cluster space as chemical symbols.")
        .def("get_cutoffs", &ClusterSpace::getCutoffs)
        .def("_get_primitive_structure", &ClusterSpace::getPrimitiveStructure)
        .def("get_multi_component_vector_permutations", &ClusterSpace::getMultiComponentVectorPermutations)
        .def("get_number_of_allowed_species_by_site", &ClusterSpace::getNumberOfAllowedSpeciesBySite)
        .def("_compute_multi_component_vectors",
             &ClusterSpace::computeMultiComponentVectors,
             "Compute the multi-component vectors (internal).")
        .def("_remove_orbits_cpp", &ClusterSpace::removeOrbits)
        .def("evaluate_cluster_function",
             &ClusterSpace::evaluateClusterFunction,
             "Evaluates value of a cluster function.")
        .def("__len__", &ClusterSpace::getClusterSpaceSize);

    py::class_<ClusterExpansionCalculator>(m, "_ClusterExpansionCalculator")
        .def(py::init<const ClusterSpace &,
                      const Structure &,
                      const double>(),
             "Initializes an icet ClusterExpansionCalculator instance.",
             py::arg("cluster_space"),
             py::arg("structure"),
             py::arg("fractional_position_tolerance"))
        .def(
            "get_cluster_vector_change",
            [](ClusterExpansionCalculator &calc,
               const py::array_t<int> &occupations,
               const int flipIndex,
               const int newOccupation)
            {
                auto cvChange = calc.getClusterVectorChange(occupations, flipIndex, newOccupation);
                return py::array(cvChange.size(), cvChange.data());
            },
            R"pbdoc(
              Returns the change in cluster vector upon flipping of one site.

              Parameters
              ----------
              occupations : list(int)
                  the occupation vector for the supercell before flip
              flip_index : int
                  local index of the supercell where flip has occured
              new_occupation : int
                  new atomic number of the flipped site
              )pbdoc",
            py::arg("occupations"),
            py::arg("flip_index"),
            py::arg("new_occupation"))
        .def(
            "get_local_cluster_vector",
            [](ClusterExpansionCalculator &calc,
               const py::array_t<int> &occupations,
               const int index)
            {
                auto localCv = calc.getLocalClusterVector(occupations, index);
                return py::array(localCv.size(), localCv.data());
            },
            R"pbdoc(
              Returns a cluster vector that only considers clusters that contain the input index.

              Parameters
              ----------
              occupations : list(int)
                  the full occupation vector for the supercell
              index : int
                  index of site whose local cluster vector should be calculated
              )pbdoc",
            py::arg("occupations"),
            py::arg("index"))
        .def(
            "get_cluster_vector",
            [](ClusterExpansionCalculator &calc,
               const py::array_t<int> &occupations)
            {
                auto cv = calc.getClusterVector(occupations);
                return py::array(cv.size(), cv.data());
            },
            R"pbdoc(
              Returns full cluster vector used in total property calculations.

              Parameters
              ----------
              occupations : list(int)
                  the occupation vector for the supercell
              )pbdoc",
            py::arg("occupations"));
}
