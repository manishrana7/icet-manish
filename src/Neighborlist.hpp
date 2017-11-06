#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <string>
#include "Structure.hpp"
#include "LatticeSite.hpp"

using namespace Eigen;

/**
class Neighborlist

Builds a (pair) neighborlist for each lattice site

*/

class Neighborlist
{
  public:
    Neighborlist(const double);
    void build(const Structure &);
    void update(const Structure &);

    std::vector<LatticeSite> getNeighbors(int index) const
    {

        if (index >= _neighbors.size() || index < 0)
        {
            throw std::out_of_range("Error: Tried accessing position at out of bound index. Neighborlist::getNeighbors");
        }
        return _neighbors[index];
    }
    ///Check if index1 and index2 are neighbors
    bool isNeighbor(const int index1, const int index2, const Vector3d offset) const
    {
        for (const auto &nbr : _neighbors[index1])
        {
            if (nbr.index() == index2) // index are the same
            {
                if (nbr.unitcellOffset() == offset) // are the offsets equal?
                {
                    return true;
                }
            }
        }
        return false;
    }

    size_t size() const
    {
        return _neighbors.size();
    }

  private:
    std::vector<std::vector<int>> latticeIndices;
    std::vector<std::vector<Vector3d>> offsets;
    std::vector<std::vector<LatticeSite>> _neighbors;
    double _cutoff;
    double DISTTOL = 1e-7;
};
