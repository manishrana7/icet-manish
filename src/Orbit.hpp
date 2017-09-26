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
    void addEquivalentSites(const std::vector<LatticeNeighbor> &latNbrs, bool sort = false)
    {
        _equivalentSites.push_back(latNbrs);
        if (sort)
        {
            sortOrbit();
        }
    }
    ///add many lattice neigbhors
    void addEquivalentSites(const std::vector<std::vector<LatticeNeighbor>> &LatticeNeighbors, bool sort = false)
    {
        _equivalentSites.insert(_equivalentSites.end(), LatticeNeighbors.begin(), LatticeNeighbors.end());
        if (sort)
        {
            sortOrbit();
        }
    }

    ///Returns amount of equivalent sites in this orbit
    size_t size() const
    {
        return _equivalentSites.size();
    }
    ///Returns the geometric size of the orbit defines as the mean distance to the center of the
    double getGeometricalSize() const
    {
        _representativeCluster.getGeometricalSize();
    }
    ///Return the sorted, reprasentative cluster for this orbit
    Cluster getRepresentativeCluster() const
    {
        return _representativeCluster;
    }

    ///Returns equivalent sites
    std::vector<std::vector<LatticeNeighbor>> getEquivalentSites() const
    {
        return _equivalentSites;
    }

    std::vector<LatticeNeighbor> GetSitesOfIndex(unsigned int index) const
    {
        if (index >= _equivalentSites.size())
        {
            throw std::out_of_range("Index out of range in function Orbit::GetSitesOfIndex");
        }
        return _equivalentSites[index];
    }
    void setEquivalentSites(const std::vector<std::vector<LatticeNeighbor>> &equivalentSites)
    {
        _equivalentSites = equivalentSites;
    }

    void sortOrbit()
    {
        std::sort(_equivalentSites.begin(), _equivalentSites.end());
    }
    ///Return the number of bodies of the cluster that represent this orbit
    unsigned int getClusterSize() const
    {
        return _representativeCluster.getNumberOfBodies();
    }
    ///Compare operator for automatic sorting in containers
    friend bool operator<(const Orbit &orbit1, const Orbit &orbit2)
    {

        // return ( orbit1.getRepresentativeCluster() < orbit2.getRepresentativeCluster());
        ///not equal size: compare by geometrical size
        if (fabs(orbit1.getGeometricalSize() - orbit2.getGeometricalSize()) > 1e-3) // @TODO: remove 1e-4 and add tunable parameter
        {
            return orbit1.getGeometricalSize() < orbit2.getGeometricalSize();
        }

        // check size of vector of equivalent sites
        if (orbit1.size() < orbit2.size())
        {
            return true;
        }
        if (orbit1.size() > orbit2.size())
        {
            return false;
        }
        //Now size of  equivalent sites vector are the same, then check the individual equivalent sites
        return orbit1.getEquivalentSites() < orbit2.getEquivalentSites();
    }
    /** 
    Returns the number of exactly equal sites in equivalent sites vector
    This is used among other things to debug orbits when duplicates is not expected
    */
    int getNumberOfDuplicates(int verbosity = 0) const;

    friend Orbit operator+(const Orbit &orbit, const Eigen::Vector3d &offset)
    {
        Orbit orbitOffset = orbit;
        for (auto &latNbrs : orbitOffset._equivalentSites)
        {
            for (auto &latNbr : latNbrs)
            {
                latNbr = latNbr + offset;
            }
        }
        return orbitOffset;
    }

  private:
    ///Representative sorted cluster for this orbit
    Cluster _representativeCluster;

    ///Container of equivalent sites for this orbit
    std::vector<std::vector<LatticeNeighbor>> _equivalentSites;
};