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

class Structure
{
  public:
    Structure(){};
    Structure(const Eigen::Matrix<double, Dynamic, 3, RowMajor> &,
              const std::vector<std::string> &,
              const Eigen::Matrix3d &,
              const std::vector<bool> &);

    double getDistance(const int, const int) const;

    /**
        Returns the distance for index1 with unitcell offset offset 1 to index2 with unit-cell offset offset2

        @TODO: use overloading here instead
    */
    double getDistance2(const int index1, const Vector3d offset1,
                        const int index2, const Vector3d offset2) const
    {
        if (index1 >= _positions.rows() or index2 >= _positions.rows())
        {
            throw std::out_of_range("Error: Tried accessing position at out of bound index. Structure::getDistance2");
        }
        Vector3d pos1 = _positions.row(index1) + offset1.transpose() * _cell;
        Vector3d pos2 = _positions.row(index2) + offset2.transpose() * _cell;
        return (pos1 - pos2).norm();
    }

    // Getters - Setters
    void setPositions(const Eigen::Matrix<double, Dynamic, 3> &positions)
    {
        _positions = positions;
    }

    ///Get position from a lattice neighbor
    Vector3d getPosition(const LatticeNeighbor &latticeNeighbor) const
    {
        if (latticeNeighbor.index >= _positions.rows() || latticeNeighbor.index < 0)
        {
            throw std::out_of_range("Error: Tried accessing position at out of bound index. Structure::getPosition");
        }
        Vector3d position = _positions.row(latticeNeighbor.index) + latticeNeighbor.unitcellOffset.transpose() * _cell;
        return position;
    }

    ///Return positions
    Eigen::Matrix<double, Dynamic, 3, RowMajor> getPositions() const
    {
        return _positions;
    }

    /**set elements using a string vector 
    @TODO: think about overloading to setElements
    */
    void setStrElements(const std::vector<std::string> &elements)
    {
        _strelements = elements;
        setElements(convertStrElements(_strelements));
    }
    ///return the elements as strings
    std::vector<std::string> getStrElements() const
    {
        return _strelements;
    }

    ///set elements using integer vector
    void setElements(const std::vector<int> &elements)
    {
        _elements = elements;
    }
    ///get elements as vector<int>
    std::vector<int> getElements() const
    {
        return _elements;
    }

    /// Get integer element of index i
    int getElement(const unsigned int i) const
    {
        if (i >= _elements.size())
        {
            std::string errorMessage = "Error: out of range in function get element:index : elements.size()  _elements.size() ";
            errorMessage += std::to_string(i) + " : ";
            errorMessage += std::to_string(_elements.size());
            throw std::out_of_range(errorMessage);
        }

        return _elements[i];
    }

    ///Set the symmetrically distinct sites using vector<int> where length of vector should match number of positions
    void setUniqueSites(const std::vector<int> &sites)
    {
        if( sites.size() != _positions.rows() )
        {
            throw std::out_of_range("Sites are not the same size as positions");            
        }
        _uniqueSites = sites;
    }
    ///Returns the symettrically distinct sites if index i and index j has the same unique site they are considered equal
    std::vector<int> getUniqueSites() const
    {
        return _uniqueSites;
    }
    ///Return the symmetrically distinct site 
    int getSite(const size_t i) const
    {
        if (i >= _uniqueSites.size())
        {
            std::string errorMessage = "Error: out of range in function getSite : index :  _uniqueSites.size() ";
            errorMessage += std::to_string(i) + " : ";
            errorMessage += std::to_string(_uniqueSites.size());

            throw std::out_of_range(errorMessage);
        }
        return _uniqueSites[i];
    }

    ///rTrue if this structure has pbc along cell axis k
    bool has_pbc(const int k) const
    {
        return _pbc[k];
    }
    ///return the pbc vector
    std::vector<bool> get_pbc() const
    {
        return _pbc;
    }
    ///set the pbc vector
    void set_pbc(const std::vector<bool> pbc)
    {
        _pbc = pbc;
    }
    ///set the lattice cell matrix
    void set_cell(const Eigen::Matrix<double, 3, 3> &cell)
    {
        _cell = cell;
    }
    ///returns the lattice cell matrix
    Eigen::Matrix<double, 3, 3> get_cell() const
    {
        return _cell;
    }

    ///Returns size of the structure i.e. number of atoms, sites, or positions in the structure
    size_t size() const
    {
        if (_elements.size() != _positions.rows())
        {
            throw std::out_of_range("Error: Positions and elements do not match in size");
        }
        return (_elements.size());
    }

    /**
    Search for the position and returns the index of the found position.

    If it does not find the position it will return -1

    argument: const double position_tolerance if the norm of difference of positions is less than this
    then equality is assumed.
    */
    int findIndexOfPosition(const Vector3d &position, const double position_tolerance = 1e-5) const
    {
        for (size_t i = 0; i < _positions.rows(); i++)
        {
            if ((_positions.row(i).transpose() - position).norm() < position_tolerance)
            {
                return i;
            }
        }

        return -1;
    }

    /**
    Finds the LatticeNeigbhor object from the position.

    The algorithm works by first extracting the fractional position.
    From the fractional position the unitcelloffset is taken by rounding the fractional coordinates to the nearest integer.
    When subtracting the fractional position with the unitcelloffset and taking the dot product with the cell 
    the remainder position is found.

    The index is found by searching for the remainder position in structure.

    if no index is found a runtime_error gets thrown.
    */
    LatticeNeighbor findLatticeNeighborFromPosition(const Vector3d &position, const double position_tolerance = 1e-5) const
    {

        ///ldlt require positive or negative semidefinite cell
        // std::cout << "position " << position << std::endl;
        //  std::cout<<"cell "<< _cell<<std::endl;
        // Vector3d position = {coordinateRound(position1[0]), coordinateRound(position1[1]), coordinateRound(position1[2]) };
        Vector3d fractional = _cell.transpose().partialPivLu().solve(position);
        // std::cout << "fractional " << fractional << std::endl;
        // Vector3d unitcellOffset = {int(round(fractional[0])), int(round(fractional[1])), int(round(fractional[2]))};
        Vector3d unitcellOffset = {int(floor(coordinateRound((double)fractional[0]))),
                                   int(floor(coordinateRound((double)fractional[1]))),
                                   int(floor(coordinateRound((double)fractional[2])))};
        // std::cout << "unitcellOffset " << unitcellOffset << std::endl;

        Vector3d remainder = (fractional - unitcellOffset).transpose() * _cell;
        // std::cout << "remainder " << remainder << std::endl;

        auto index = findIndexOfPosition(remainder, position_tolerance);
        if (index == -1)
        {
            std::cout<<"Positions"<<std::endl;
           for(int i =0; i < size(); i++)                
           {
               Vector3d pos = _positions.row(i);
               std::cout<<pos[0]<<" "<<pos[1]<< " "<<pos[2]<<std::endl;
           }
           std::cout<<"Positions done"<<std::endl;
        ///ldlt require positive or negative semidefinite cell
        std::cout << "position " << position << std::endl;
        std::cout<<"cell "<< _cell<<std::endl;
        std::cout << "fractional " << fractional << std::endl;
        Vector3d unitcellOffset = {int(round(fractional[0])), int(round(fractional[1])), int(round(fractional[2]))};
        std::cout << "unitcellOffset " << unitcellOffset << std::endl;
        std::cout << "remainder " << remainder << std::endl;

            throw std::runtime_error("Did not find position in function findLatticeNeighborFromPosition in Structure");
        }

        LatticeNeighbor ret = LatticeNeighbor(index, unitcellOffset);
        return ret;
    }

    ///Round to nearest integer toward zero
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
    /**
    Finds a vector of lattice neigbhors from a vector of positions

    */
    std::vector<LatticeNeighbor> findLatticeNeighborsFromPositions(const std::vector<Vector3d> &positions, const double position_tolerance = 1e-5) const
    {
        std::vector<LatticeNeighbor> latNbrVector;
        latNbrVector.reserve(positions.size());

        for (const Vector3d position : positions)
        {
            latNbrVector.push_back(findLatticeNeighborFromPosition(position, position_tolerance));
        }

        return latNbrVector;
    }
    
    /// Return number of allowed components on site i
    int getMi(const unsigned int i) const
    {
        if (i >= _allowedComponents.size())
        {
            std::string errorMessage = "Error: out of range in function getMi : index :  _allowedComponents.size() ";
            errorMessage += std::to_string(i) + " : ";
            errorMessage += std::to_string(_allowedComponents.size());

            throw std::out_of_range(errorMessage);
        }
        return _allowedComponents[i];
    }
    
    ///Set allowed components on each site
    void setAllowedComponents(const std::vector<int> &allowedComponents)
    {
        _allowedComponents = allowedComponents;
    }   

  private:
    Eigen::Matrix<double, Dynamic, 3, RowMajor> _positions;
    Eigen::Matrix3d _cell;
    std::vector<int> _elements;
    std::vector<std::string> _strelements;
    std::vector<bool> _pbc;
    std::vector<int> _uniqueSites;
    std::vector<int> _allowedComponents;
    
    std::vector<int> convertStrElements(const std::vector<std::string> &elements)
    {
        std::vector<int> intElements(elements.size());
        for (int i = 0; i < elements.size(); i++)
        {
            intElements[i] = PeriodicTable::strInt[elements[i]];
        }
        return intElements;
    }

    std::vector<std::string> convertIntElements(const std::vector<int> &elements)
    {
        std::vector<std::string> strElements(elements.size());
        for (int i = 0; i < elements.size(); i++)
        {
            strElements[i] = PeriodicTable::intStr[elements[i]];
        }
        return strElements;
    }
    ///rounds val to precision
    double coordinateRound(const double &val) const
    {
        double precision = 1e-7;
        return round(val * 1.0 / precision) / (1.0 / precision);
    }
};
