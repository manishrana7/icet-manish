#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "PeriodicTable.hpp"
#include "LatticeNeighbor.hpp"

using namespace Eigen;
namespace py = pybind11;

/**
  @brief Class for storing a structure.
  @details This class stores the cell metric, positions, chemical symbols, and
  periodic boundary conditions. It also holds information pertaining to the
  components that are allowed on each site. It also provides functionality for
  computing e.g., distances between sites.
  @todo rename element/elements/_elements to atomicNumber/atomicNumbers/_atomicNumbers
*/
class Structure
{
  public:

    /// Constructor.
    Structure(){};

    /// Constructor.
    Structure(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &,
              const std::vector<std::string> &,
              const Eigen::Matrix3d &,
              const std::vector<bool> &,
              double);

    /// Return distance vector between two sites.
    double getDistance(const int, const int) const;

    /// Return distance vector between two sites.
    double getDistance2(const int, const Vector3d, const int, const Vector3d) const;

    /// Return the position of a site in Cartesian coordinates.
    Vector3d getPosition(const LatticeNeighbor &) const;

    /// Return atomic number of site.
    int getElement(const unsigned int) const;

    /// Return the list of unique sites.
    std::vector<int> getUniqueSites() const { return _uniqueSites; }

    /// Set list of unique sites.
    /// @todo add example for how the unique sites are supposed to work.
    void setUniqueSites(const std::vector<int> &);

    /// Return a unique site.
    /// @todo rename to getUniqueSite
    int getSite(const size_t) const;

    /// Return index of site that matches the given position.
    /// @todo rename to findSiteByPosition
    int findIndexOfPosition(const Vector3d &) const;

    /// Return LatticeNeighbor object that matches the given position.
    /// @todo rename to findLatticeNeighborByPosition
    LatticeNeighbor findLatticeNeighborFromPosition(const Vector3d &) const;

    /// Return list of LatticeNeighbor objects that matche a given list of positions.
    std::vector<LatticeNeighbor> findLatticeNeighborsFromPositions(const std::vector<Vector3d> &) const;

    /// Representation of class as string (cast).
    operator std::string () const
    {
        std::string str;
        str += "Structure";
        str += " nsites: " + std::to_string(_positions.size());
        return str;
    }

  public:

    /// Return the size of the structure, i.e., the number of sites.
    size_t size() const { return (_elements.size()); }

    // Set the atomic positions.
    void setPositions(const Eigen::Matrix<double, Dynamic, 3> &positions) { _positions = positions; }

    /// Return positions.
    Eigen::Matrix<double, Dynamic, 3, RowMajor> getPositions() const { return _positions; }

    /// Set elements.
    /// @todo Rename to setAtomicNumbers
    /// @todo Rename to setChemicalSymbols
    /// @todo Think about overloading to setElements
    void setElements(const std::vector<int> &elements) { _elements = elements; }
    void setStrElements(const std::vector<std::string> &elements)
    {
        _strelements = elements;
        setElements(convertStrElements(_strelements));
    }

    /// Return elements.
    std::vector<int> getElements() const { return _elements; }
    std::vector<std::string> getStrElements() const { return _strelements; }

    /**
      @brief Return periodic boundary condition along a given direction.
      @details Return True if periodic boundary conditions are active along direction.
      @param k index to direction [0, 1, 2]
    **/
    bool has_pbc(const int k) const { return _pbc[k]; }

    /// Return periodic boundary conditions.
    std::vector<bool> get_pbc() const { return _pbc; }

    /// Set periodic boundary conditions.
    void set_pbc(const std::vector<bool> pbc) { _pbc = pbc; }

    /// Set the cell metric.
    /// @todo switch to CamelCase
    void set_cell(const Eigen::Matrix<double, 3, 3> &cell) { _cell = cell; }

    /// Return the cell metric.
    /// @todo switch to CamelCase
    Eigen::Matrix<double, 3, 3> get_cell() const { return _cell; }

    /// Set allowed components for each site.
    /// @todo rename to setNumberOfAllowedComponents
    void setAllowedComponents(const std::vector<int> &);
    void setAllowedComponents(const int);

    /// Return number of allowed components on each site.
    /// @todo rename to getNumberOfAllowedComponents
    int getMi(const unsigned int) const;

    /// Set tolerance.
    void setTolerance(double tolerance) { _tolerance = tolerance; }

    /// Return tolerance.
    double getTolerance() const { return _tolerance; }


  private:

    /**
      @brief Convert chemical symbols to atomic numbers.
      @param elements vector of strings to be converted
      @todo rename to convertChemicalSymbolsToAtomicNumbers
    **/
    std::vector<int> convertStrElements(const std::vector<std::string> &elements)
    {
        std::vector<int> intElements(elements.size());
        for (int i = 0; i < elements.size(); i++)
        {
            intElements[i] = PeriodicTable::strInt[elements[i]];
        }
        return intElements;
    }

    /**
      @brief Convert chemical symbols to atomic numbers.
      @param elements vector of strings to be converted
      @todo rename to convertAtomicNumbersToChemicalSymbols
    **/
    std::vector<std::string> convertIntElements(const std::vector<int> &elements)
    {
        std::vector<std::string> strElements(elements.size());
        for (int i = 0; i < elements.size(); i++)
        {
            strElements[i] = PeriodicTable::intStr[elements[i]];
        }
        return strElements;
    }

    /// Round float number to given tolerance.
    /// @todo rename to roundFloat
    /// @todo move to a more general location.
    double coordinateRound(const double &val, const double rounding_tolerance = 1e-7) const
    {
        return round(val * 1.0 / rounding_tolerance) / (1.0 / rounding_tolerance);
    }

    /// Round to nearest integer toward zero.
    /// @todo move to a more general location.
    int nearestIntegerTowardZero(const double value) const
    {
        if (value > 0)
        {
            return int(floor(value));
        }
        else
        {
            return int(floor(value));
        }
    }

  private:

    /// Positions of sites.
    Eigen::Matrix<double, Dynamic, 3, RowMajor> _positions;

    /// Cell metric.
    Eigen::Matrix3d _cell;

    /// List of atomic numbers.
    std::vector<int> _elements;

    /// List of chemical symbols
    /// @todo get rid of it; the information should only be generated on-demand
    std::vector<std::string> _strelements;

    /// Periodic boundary conditions.
    std::vector<bool> _pbc;

    /// List of unique sites.
    /// @todo currently not used
    std::vector<int> _uniqueSites;

    /// List of the number of allowed components on each site.
    std::vector<int> _numbersOfAllowedComponents;

    /// tolerance used for rounding positions.
    double _tolerance;

};
