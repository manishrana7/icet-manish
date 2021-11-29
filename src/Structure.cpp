#include "Structure.hpp"

using namespace Eigen;

/**
  @details Initializes an icet Structure instance.
  @param positions list of positions in Cartesian coordinates
  @param atomicNumbers list of atomic numbers
  @param cell cell metric
  @param pbc periodic boundary conditions
**/
Structure::Structure(const Matrix<double, Dynamic, 3, RowMajor> &positions,
                     const py::array_t<int> &atomicNumbers,
                     const Matrix3d &cell,
                     const std::vector<bool> &pbc)
    : _atomicNumbers(atomicNumbers), _cell(cell), _pbc(pbc), _positions(positions)
{
}

/**
  @details This function returns the position of a site.
  @param latticeNeighbor site for which to obtain the position
  @returns a 3-dimensional position vector
*/
Vector3d Structure::position(const LatticeSite &latticeNeighbor) const
{
    if (latticeNeighbor.index() >= (size_t)_positions.rows())
    {
        std::ostringstream msg;
        msg << "Site index out of bounds";
        msg << " index: " << latticeNeighbor.index();
        msg << " number of positions: " << _positions.rows();
        msg << " (Structure::position)";
        throw std::out_of_range(msg.str());
    }
    Vector3d position = _positions.row(latticeNeighbor.index()) + latticeNeighbor.unitcellOffset().transpose().cast<double>() * _cell;
    return position;
}
/**
@details This function returns the position of a specific site in Cartesian coordinates.
@param index index of the site
 **/
Vector3d Structure::positionByIndex(const size_t &index) const
{
    Vector3d position = _positions.row(index);
    return position;
}
/**
  @details This function returns the atomic number of a site.
  @param index index of site
  @returns atomic number
**/
int Structure::getAtomicNumber(const size_t index) const
{
    if (index >= _atomicNumbers.size())
    {
        std::ostringstream msg;
        msg << "Site index out of bounds";
        msg << " index: " << index;
        msg << " nsites: " << _atomicNumbers.size();
        msg << " (Structure::getAtomicNumber)";
        throw std::out_of_range(msg.str());
    }
    return _atomicNumbers.at(index);
}

/**
  @details This function returns the LatticeSite object the position of
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

  @param position position to match in Cartesian coordinates
  @param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates

  @returns LatticeSite object
*/
LatticeSite Structure::findLatticeSiteByPosition(const Vector3d &position, const double fractionalPositionTolerance) const
{
    /// Loop over all positions
    for (size_t i = 0; i < (size_t)_positions.rows(); i++)
    {
        Vector3d distanceVector = position - _positions.row(i).transpose();
        Vector3d fractionalDistanceVector = _cell.transpose().partialPivLu().solve(distanceVector);

        // Check whether whether the fractionalDistanceVector is all integers,
        // if it is we have found the corresponding lattice site
        Vector3d latticeVector = {round(fractionalDistanceVector[0]),
                                  round(fractionalDistanceVector[1]),
                                  round(fractionalDistanceVector[2])};
        if ((fractionalDistanceVector - latticeVector).norm() < fractionalPositionTolerance)
        {
            return LatticeSite(i, latticeVector.cast<int>());
        }
    }

    Vector3d fractionalPosition = _cell.transpose().partialPivLu().solve(position);
    std::ostringstream msg;
    msg << "Failed to find site by position (findLatticeSiteByPosition)." << std::endl;
    msg << "Try increasing symprec or position_tolerance." << std::endl;
    msg << "position: " << position[0] << " " << position[1] << " " << position[2] << std::endl;
    msg << "fractional position: " << fractionalPosition[0] << " " << fractionalPosition[1] << " " << fractionalPosition[2] << std::endl;
    msg << "fractional position tolerance: " << fractionalPositionTolerance;
    throw std::runtime_error(msg.str());
}

/**
  @brief Prescribes the atomic numbers that are allowed on each site in the structure.
  @param atomicNumbers
    Nested list with the atomic numbers allowed on each site in the structure.
**/
void Structure::setAllowedAtomicNumbers(const std::vector<std::vector<int>> &atomicNumbers)
{
    if (atomicNumbers.size() != size())
    {
        std::ostringstream msg;
        msg << "Size of input list incompatible with structure";
        msg << " length: " << atomicNumbers.size();
        msg << " nsites: " << size();
        msg << " (Structure::setAllowedAtomicNumbers)";
        throw std::out_of_range(msg.str());
    }
    _allowedAtomicNumbers = atomicNumbers;
}

/**
  @details This function returns the number of components allowed on a
  given site.
  @param index index of the site
  @returns the number of the allowed components
**/
int Structure::getNumberOfAllowedSpeciesBySite(const size_t index) const
{
    if (!hasAllowedAtomicNumbers())
    {
        std::ostringstream msg;
        msg << "Allowed atomic numbers per site not set in this structure";
        msg << " (Structure::getNumberOfAllowedSpeciesBySite)";
        throw std::out_of_range(msg.str());
    }
    else if (index >= _allowedAtomicNumbers.size())
    {
        std::ostringstream msg;
        msg << "Site index out of bounds";
        msg << " index: " << index;
        msg << " nsites: " << _allowedAtomicNumbers.size();
        msg << " (Structure::getNumberOfAllowedSpeciesBySite)";
        throw std::out_of_range(msg.str());
    }
    return _allowedAtomicNumbers[index].size();
}

/**
  @details This function returns the a vector with number of components allowed on each site index
  @param sites indices of sites
  @returns the list of number of allowed components for each site
**/
std::vector<int> Structure::getNumberOfAllowedSpeciesBySites(const std::vector<LatticeSite> &sites) const
{
    std::vector<int> numberOfAllowedSpecies(sites.size());
    int i = -1;
    for (const auto site : sites)
    {
        i++;
        numberOfAllowedSpecies[i] = getNumberOfAllowedSpeciesBySite(site.index());
    }
    return numberOfAllowedSpecies;
}
