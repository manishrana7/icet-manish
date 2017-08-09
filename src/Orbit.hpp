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

    void addEquivalentSites(const std::vector<LatticeNeighbor> &latNbrs)
    {
        _equivalentSites.push_back(latNbrs);
    }



    friend bool operator<(const Orbit &orbit1, const Orbit &orbit2)
    {
        if( orbit1.getRepresentativeCluster() < orbit2.getRepresentativeCluster())
        {
            return true;
        }
        if( orbit1.getRepresentativeCluster() > orbit2.getRepresentativeCluster())
        {
            return false;
        }
        //representative cluster is equal
        //Try comparing length of equivalent sites
        if( orbit1.size() < orbit2.size() )
        {
            return true;
        }
        if( orbit1.size() > orbit2.size() )
        {
            return false;
        }   
        //Both representative cluster and size of equivalent sites are equal.
        //throw error to see if this ever happens

        throw std::runtime_error("Both representative cluster and size of equivalent sites are equal in orbit < comparison");
    }

    size_t size() const
    {
        return _equivalentSites.size();
    }

    Cluster getRepresentativeCluster() const
    {
        return _sortedCluster;
    }
    std::vector<<std::vector<LatticeNeighbor>> getEquivalentSites() const
    {
        return _equivalentSites;
    }
    private:
    Cluster _sortedCluster;
    std::vector<std::vector<LatticeNeighbor>> _equivalentSites;
};