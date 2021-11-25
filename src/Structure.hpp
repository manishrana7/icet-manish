#pragma once

#include <Eigen/Dense>
#include <pybind11/numpy.h>

#include "LatticeSite.hpp"

using namespace Eigen;
namespace py = pybind11;

/**
  @brief Class for storing a structure.
  @details This class stores the cell metric, positions, atomic numbers, and
  periodic boundary conditions that describe a structure. It also holds
  information pertaining to the components that are allowed on each site and
  provides functionality for computing distances between sites.
*/
class Structure
{
public:
  /// Default constructor.
  Structure(){};

  /// Overloaded constructor.
  Structure(const Matrix<double, Dynamic, 3, RowMajor> &,
            const py::array_t<int> &,
            const Matrix3d &,
            const std::vector<bool> &);

  /// Returns distance vector between two sites.
  double getDistance(const size_t, const size_t, const Vector3i, const Vector3i) const;

  /// Return the position of a site in Cartesian coordinates.
  Vector3d getPosition(const LatticeSite &) const;

  /// Return the position of a site in Cartesian coordinates.
  Vector3d getPositionByIndex(const size_t &index) const;

  /// Returns atomic number of site.
  int getAtomicNumber(const size_t) const;

  /// Return LatticeSite object that matches the given position.
  LatticeSite findLatticeSiteByPosition(const Vector3d &, const double) const;

  /// Returns the size of the structure, i.e., the number of sites.
  size_t size() const { return _atomicNumbers.size(); }

  /// Set the atomic positions.
  void setPositions(const Matrix<double, Dynamic, 3> &positions) { _positions = positions; }

  /// Returns positions.
  Matrix<double, Dynamic, 3, RowMajor> getPositions() const { return _positions; }

  /// Set atomic numbers.
  void setAtomicNumbers(const py::array_t<int> &atomicNumbers) { _atomicNumbers = atomicNumbers; }

  /// Returns atomic numbers.
  const py::array_t<int> &getAtomicNumbers() const { return _atomicNumbers; }

  /// Returns periodic boundary conditions.
  std::vector<bool> getPBC() const { return _pbc; }

  /// Set periodic boundary conditions.
  void setPBC(const std::vector<bool> pbc) { _pbc = pbc; }

  /// Set the cell metric.
  void setCell(const Matrix<double, 3, 3> &cell) { _cell = cell; }

  /// Returns the cell metric.
  Matrix<double, 3, 3> getCell() const { return _cell; }

  /// Set allowed components for each site by vector.
  void setNumberOfAllowedSpecies(const std::vector<int> &);

  /// Returns number of allowed components on each site.
  int getNumberOfAllowedSpeciesBySite(const size_t) const;

  /// Returns number of allowed components on each site.
  std::vector<int> getNumberOfAllowedSpeciesBySites(const std::vector<LatticeSite> &) const;

private:
  /// List of atomic numbers.
  py::array_t<int> _atomicNumbers;

  /// Positions of sites in Cartesian coordinates.
  Matrix<double, Dynamic, 3, RowMajor> _positions;

  /// Cell metric.
  Matrix3d _cell;

  /// Periodic boundary conditions.
  std::vector<bool> _pbc;

  /// List of the number of allowed components on each site.
  std::vector<int> _numbersOfAllowedSpecies;
};
