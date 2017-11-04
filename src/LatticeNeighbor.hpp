#pragma once
#include <eigen3/Eigen/Dense>
#include <boost/functional/hash.hpp>
#include <iostream>
using boost::hash;
using boost::hash_combine;
using boost::hash_value;

/**
@brief Class for storing information concerning a lattice site.
@details This class provides functionality for handling
*/
struct LatticeNeighbor
{

public:
    /// Empty constructor.
    LatticeNeighbor() { }

    /**
    @brief Constructor.
    @param index site index
    @param unitcellOffset offset of site relative to unit cell at origin in units of lattice vectors
    */
    LatticeNeighbor(const int index, const Eigen::Vector3d unitcellOffset)
    {
        _index = index;
        _unitcellOffset = unitcellOffset;
    }

    /// Return index of site.
    int index() const { return _index; }

    /// Set index of site.
    void setIndex(int index) { _index = index; }

    /// Return offset relative to unit cell at origin in units of lattice vectors.
    Eigen::Vector3d unitcellOffset() const { return _unitcellOffset; }

    /// Set offset relative to unit cell at origin in units of lattice vectors.
    void setUnitcellOffset(Eigen::Vector3d offset) { _unitcellOffset = offset; }

    /// Add offset relative to unit cell at origin in units of lattice vectors.
    void addUnitcellOffset(Eigen::Vector3d offset) { _unitcellOffset += offset; }

    /// Smaller than operator.
    bool operator<(const LatticeNeighbor &other) const
    {
        if (_index == other.index())
        {
            for (int i = 0; i < 3; i++)
            {
                if (_unitcellOffset[i] != other.unitcellOffset()[i])
                {
                    return _unitcellOffset[i] < other.unitcellOffset()[i];
                }
            }
        }
        return _index < other.index();
    }

    /// Equality operator.
    bool operator==(const LatticeNeighbor &other) const
    {
        if (_index != other.index())
        {
            return false;
        }

        for (int i = 0; i < 3; i++)
        {
            if (_unitcellOffset[i] != other.unitcellOffset()[i])
            {
                return false;
            }
        }
        return true;
    }

    /// Addition operator.
    friend LatticeNeighbor operator+(const LatticeNeighbor &latticeNeighbor, const Eigen::Vector3d &offset)
    {
        LatticeNeighbor latnbr = LatticeNeighbor(latticeNeighbor.index(), latticeNeighbor.unitcellOffset() + offset);
        return latnbr;
    }

    /// Substraction operator.
    friend LatticeNeighbor operator-(const LatticeNeighbor &latticeNeighbor, const Eigen::Vector3d &offset)
    {
        LatticeNeighbor latnbr = LatticeNeighbor(latticeNeighbor.index(), latticeNeighbor.unitcellOffset() - offset);
        return latnbr;
    }

    /// Addition and assignment operator.
    friend LatticeNeighbor operator+=(const LatticeNeighbor &latticeNeighbor, const Eigen::Vector3d &offset)
    {
        LatticeNeighbor latnbr = LatticeNeighbor(latticeNeighbor.index(), latticeNeighbor.unitcellOffset() + offset);
        return latnbr;
    }

    /// Cast object as string.
    operator std::string () const
    {
        std::string str = std::to_string(_index) + " :";
        for (int i = 0; i < 3; i++)
        {
            str += " " + std::to_string(_unitcellOffset[i]);
        }
        return str;
    }

    /// Write class information to stdout.
    void print() const
    {
        std::cout << _index << " : ";
        for (int i = 0; i < 3; i++)
        {
            std::cout << _unitcellOffset[i] << " ";
        }
        std::cout << std::endl;
    }

private:

    /// Site index.
    int _index;

    /// Offset relative to the unit cell at the origin in units of lattice vectors.
    Eigen::Vector3d _unitcellOffset;

};

namespace std
{
    template <>
    struct hash<LatticeNeighbor>
    {
        size_t
        operator()(const LatticeNeighbor &k) const
        {

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            size_t seed = 0;
            hash_combine(seed, hash_value(k.index()));

            for (int i = 0; i < 3; i++)
            {
                hash_combine(seed, hash_value(k.unitcellOffset()[i]));
            }
            return seed;
        }
    };

    template <>
    struct hash<std::vector<LatticeNeighbor>>
    {
        size_t
        operator()(const std::vector<LatticeNeighbor> &k) const
        {

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            size_t seed = 0;
            for (const auto &latNbr : k)
            {
                hash_combine(seed, std::hash<LatticeNeighbor>{}(latNbr));
            }
            return seed;
        }
    };
}
