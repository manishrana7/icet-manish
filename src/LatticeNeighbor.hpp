#pragma once
#include <eigen3/Eigen/Dense>
#include <boost/functional/hash.hpp>
#include <iostream>
using boost::hash;
using boost::hash_combine;
using boost::hash_value;

struct LatticeNeighbor
{

    LatticeNeighbor()
    {
        //empty constructor
    }

    LatticeNeighbor(const int newIndex, const Eigen::Vector3d newUnitcellOffset)
    {
        index = newIndex;
        unitcellOffset = newUnitcellOffset;
    }   

    int index;
    Eigen::Vector3d unitcellOffset;
    bool operator<(const LatticeNeighbor &other) const
    {
        if (index == other.index)
        {
            for (int i = 0; i < 3; i++)
            {
                if (unitcellOffset[i] != other.unitcellOffset[i])
                {
                    return unitcellOffset[i] < other.unitcellOffset[i];
                }
            }
        }
        return index < other.index;
    }
    bool operator==(const LatticeNeighbor &other) const
    {
        if (index != other.index)
        {
            return false;
        }

        for (int i = 0; i < 3; i++)
        {
            if (unitcellOffset[i] != other.unitcellOffset[i])
            {
                return false;
            }
        }
        return true;
    }

    void print() const
    {
        std::cout << index << " : ";
        for (int i = 0; i < 3; i++)
        {
            std::cout << unitcellOffset[i] << " ";
        }
        std::cout << std::endl;
    }
};

#include <boost/functional/hash.hpp>

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
        hash_combine(seed, hash_value(k.index));

        for (int i = 0; i < 3; i++)
        {
            hash_combine(seed, hash_value(k.unitcellOffset[i]));
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
