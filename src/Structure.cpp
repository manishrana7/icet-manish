#include "Structure.hpp"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace Eigen;
namespace py = pybind11;


/// Constructor.
Structure::Structure(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &positions,
                     const std::vector<std::string> &elements,
                     const Eigen::Matrix3d &cell,
                     const std::vector<bool> &pbc,
                     double tolerance)
{
    setPositions(positions);
    setStrElements(elements);
    _cell = cell;
    _pbc = pbc;
    _tolerance = tolerance;
    _uniqueSites.resize(elements.size());
    _numbersOfAllowedComponents.resize(positions.rows());
}

/**
  @details This function computes the distance between two sites.
  @param index1 index of the first site
  @param index2 index of the second site
  @todo Use overloading here instead
  @todo Add  mic functionality
*/
double Structure::getDistance(const int index1, const int index2) const
{

    if (index1 >= _positions.rows() or index2 >= _positions.rows())
    {
        throw std::out_of_range("Error: Tried accessing position at out of bound index. Structure::getDistance");
    }

    const Eigen::Vector3d pos1 = _positions.row(index1);
    const Eigen::Vector3d pos2 = _positions.row(index2);

    return (pos1 - pos2).norm();
}

/**
  @details This function computes the distance between two sites.
  @param index1 index of the first site
  @param index2 index of the second site
  @param offset1 offset of site 1 relative to origin in units of lattice vectors
  @param offset2 offset of site 2 relative to origin in units of lattice vectors
  @todo Use overloading here instead
*/
double Structure::getDistance2(const int index1, const Vector3d offset1,
                    const int index2, const Vector3d offset2) const
{
    if (index1 < 0 || index1 >= _positions.rows() or
        index2 < 0 || index2 >= _positions.rows())
    {
        std::string errorMessage = "Site index out of bounds ";
        errorMessage += " index1:" + std::to_string(index1);
        errorMessage += " index2:" + std::to_string(index1);
        errorMessage += " npositions: " + std::to_string(_positions.rows());
        errorMessage += " (Structure::getDistance2)";
        throw std::out_of_range(errorMessage);
    }
    Vector3d pos1 = _positions.row(index1) + offset1.transpose() * _cell;
    Vector3d pos2 = _positions.row(index2) + offset2.transpose() * _cell;
    return (pos1 - pos2).norm();
}

/**
  @details This function returns the position of a site.
  @param latticeNeighbor site for which to obtain the position
  @returns a 3-dimensional position vector
  @todo rename LatticeNeighbor class
*/
Vector3d Structure::getPosition(const LatticeNeighbor &latticeNeighbor) const
{
    if (latticeNeighbor.index < 0 || latticeNeighbor.index >= _positions.rows())
    {
        std::string errorMessage = "Site index out of bounds";
        errorMessage += " index: " + std::to_string(latticeNeighbor.index);
        errorMessage += " npositions: " + std::to_string(_positions.rows());
        errorMessage += " (Structure::getPosition)";
        throw std::out_of_range(errorMessage);
    }
    Vector3d position = _positions.row(latticeNeighbor.index) + latticeNeighbor.unitcellOffset.transpose() * _cell;
    return position;
}

/**
  @details This function returns the atomic number of a site.
  @param i index of site
  @returns atomic number
  @todo rename to getAtomicNumber
**/
int Structure::getElement(const unsigned int i) const
{
    if (i >= _elements.size())
    {
        std::string errorMessage = "Site index out of bounds";
        errorMessage += " i: " + std::to_string(i);
        errorMessage += " nsites: " + std::to_string(_elements.size());
        errorMessage += " (Structure::getElement)";
        throw std::out_of_range(errorMessage);
    }
    return _elements[i];
}

/**
  @details This function sets the symmetrically distinct sites associated
  with the structure. It requires a vector<int> as input the length of
  which  must match the number of positions.
  @param sites list of integers
**/
void Structure::setUniqueSites(const std::vector<int> &sites)
{
    if (sites.size() != _positions.rows())
    {
        std::string errorMessage = "Length of input vector does not match number of sites";
        errorMessage += " nsites: " + std::to_string(sites.size());
        errorMessage += " npositions: " + std::to_string(_positions.rows());
        errorMessage += " (Structure::setUniqueSites)";
        throw std::out_of_range(errorMessage);
    }
    _uniqueSites = sites;
}

/**
  @details This function returns the index of a unique site from the list of
  unique sites.
  @param i index of site
  @returns index of unique site
**/
int Structure::getSite(const size_t i) const
{
    if (i < 0 || i >= _uniqueSites.size())
    {
        std::string errorMessage = "Site index out of bounds";
        errorMessage += " i: " + std::to_string(i);
        errorMessage += " nsites: " + std::to_string(_uniqueSites.size());
        errorMessage += " (Structure::getSite)";
        throw std::out_of_range(errorMessage);
    }
    return _uniqueSites[i];
}

/**
  @details This function returns the index of the site the position of
  which matches the input position to the tolerance specified for this
  structure. If the function failes to find a match it returns -1

  @param position position vector to match

  @returns index of site
**/
int Structure::findIndexOfPosition(const Vector3d &position) const
{
    for (size_t i = 0; i < _positions.rows(); i++)
    {
        if ((_positions.row(i).transpose() - position).norm() < _tolerance)
        {
            return i;
        }
    }

    return -1;
}

/**
  @details This function returns the LatticeNeighbor object the position of
  which matches the input position to the tolerance specified for this
  structure.

  The algorithm commences by extracting the fractional position.
  From the fractional position the unitcell offset is taken by rounding the
  fractional coordinates to the nearest integer.
  By subtracting the unitcell offset from the fractional position and taking
  the dot product with the cell the position relative to the primitive cell is
  found.
  The index is found by searching for the remainder position in structure.
  If no index is found a runtime_error is thrown.

  @param position position vector to match

  @returns LatticeNeighbor object
*/
LatticeNeighbor Structure::findLatticeNeighborFromPosition(const Vector3d &position) const
{

    Vector3d fractional = _cell.transpose().partialPivLu().solve(position);
    Vector3d unitcellOffset;
    for (int i = 0; i < 3 ; i++)
    {
        unitcellOffset[i] = floor(coordinateRound((double)fractional[i]));
        if ( fabs(unitcellOffset[i] - fractional[i]) > 1.0 - _tolerance && has_pbc(i))
        {
            unitcellOffset[i] = int(round(fractional[i]));
        }
    }
    Vector3d remainder = (fractional - unitcellOffset).transpose() * _cell;

    auto index = findIndexOfPosition(remainder);
    if (index == -1)
    {
        std::string errorMessage = "Failed to find site by position (findLatticeNeighborFromPosition)";
        throw std::runtime_error(errorMessage);
    }

    LatticeNeighbor ret = LatticeNeighbor(index, unitcellOffset);
    return ret;
}

/**
  @details This function returns a list ofLatticeNeighbor object the position
  of each matches the respective entry in the list of input positions to the
  tolerance specified for this structure. Internally this function uses
  Structure::findLatticeNeighborFromPosition.

  @param positions list of position vectors to match

  @returns list of LatticeNeighbor objects
*/
std::vector<LatticeNeighbor> Structure::findLatticeNeighborsFromPositions(const std::vector<Vector3d> &positions) const
{
    std::vector<LatticeNeighbor> latNbrVector;
    latNbrVector.reserve(positions.size());

    for (const Vector3d position : positions)
    {
        latNbrVector.push_back(findLatticeNeighborFromPosition(position));
    }

    return latNbrVector;
}

/**
  @details This function allows one to specify the number of components
  that are allowed on each lattice site. This can be employed to construct
  "parallel" cluster expansions such as in (A,B) on site #1 with h(C,D) on
  site #2.
  @param numbersOfAllowedComponents list with the number of components
  allowed on each site
**/
void Structure::setAllowedComponents(const std::vector<int> &numbersOfAllowedComponents)
{
    if (numbersOfAllowedComponents.size() != size())
    {
        std::string errorMessage;
        errorMessage += "Size of input list incompatible with structure";
        errorMessage += " length: " + std::to_string(numbersOfAllowedComponents.size());
        errorMessage += " nsites: " + std::to_string(size());
        errorMessage += " (Structure::setAllowedComponents)";
        throw std::out_of_range(errorMessage);
    }
    _numbersOfAllowedComponents = numbersOfAllowedComponents;
}

/**
  @details This function allows one to specify the number of components
  that are allowed on each lattice site. This can be employed to construct
  "parallel" cluster expansions such as in (A,B) on site #1 with h(C,D) on
  site #2.
  @param numberOfAllowedComponents number of components allowed
**/
void Structure::setAllowedComponents(const int numberOfAllowedComponents)
{
    std::vector<int> numbersOfAllowedComponents(_elements.size(), numberOfAllowedComponents);
    _numbersOfAllowedComponents = numbersOfAllowedComponents;
}

/**
  @details This function returns the number of components allowed on a
  given site.
  @param i index of the site
  @returns the number of the allowed components
**/
int Structure::getMi(const unsigned int i) const
{
    if (i >= _numbersOfAllowedComponents.size())
    {
        std::string errorMessage = "Site index out of bounds";
        errorMessage += " i: " + std::to_string(i);
        errorMessage += " nsites: " + std::to_string(_numbersOfAllowedComponents.size());
        errorMessage += " (Structure::getMi)";
        throw std::out_of_range(errorMessage);
    }
    return _numbersOfAllowedComponents[i];
}
