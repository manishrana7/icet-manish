#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <vector>
#include <string>
#include "LatticeNeighbor.hpp"
#include "Cluster.hpp"
using namespace Eigen;

/**
Class Orbit

contains equivalent vector<LatticeNeighbors>
contains a sorted Cluster for representation

Can be compared to other orbits

*/

class Orbit
{
  public:
    Orbit(const Cluster &);

    ///Add a group sites that are equivalent to the ones in this orbit
    void addEquivalentSites(const std::vector<LatticeNeighbor> &latNbrs)
    {
        _equivalentSites.push_back(latNbrs);
    }


    ///Returns amount of equivalent sites in this orbit
    size_t size() const
    {
        return _equivalentSites.size();
    }

    ///Return the sorted, reprasentative cluster for this orbit
    Cluster getRepresentativeCluster() const
    {
        return _representativeCluster;
    }

    ///Returns equivalent sites
    std::vector<std::vector<LatticeNeighbor>>  getEquivalentSites() const
    {
        return _equivalentSites;
    }


      ///Compare operator for automatic sorting in containers
    friend bool operator<(const Orbit &orbit1, const Orbit &orbit2)
    {
        if (orbit1.getRepresentativeCluster() < orbit2.getRepresentativeCluster())
        {
            return true;
        }
        //note the order is changed here "o2 < o1"
        if (orbit2.getRepresentativeCluster() < orbit1.getRepresentativeCluster())
        {
            return false;
        }
        //representative cluster is equal
        //Try comparing length of equivalent sites
        if (orbit1.size() < orbit2.size())
        {
            return true;
        }
        if (orbit1.size() > orbit2.size())
        {
            return false;
        }
        //Both representative cluster and size of equivalent sites are equal.
        //throw error to see if this ever happens

        throw std::runtime_error("Both representative cluster and size of equivalent sites are equal in orbit < comparison");
    }

    int getNumberOfDuplicates(int verbosity=0) const;
  private:
    ///Reprasentative sorted cluster for this orbit
    Cluster _representativeCluster;

    ///Container of equivalent sites for this orbit
    std::vector<std::vector<LatticeNeighbor>> _equivalentSites;
};