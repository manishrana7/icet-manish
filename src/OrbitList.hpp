#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <vector>
#include "Orbit.hpp"
#include "ManybodyNeighborlist.hpp"
#include "Structure.hpp"
#include "Cluster.hpp"
#include <unordered_map>
/**
Class OrbitList

contains a sorted vector or orbits


*/

class OrbitList
{
  public:
    OrbitList();
    OrbitList(const ManybodyNeighborlist &, const Structure &);

    ///Add a group sites that are equivalent to the ones in this orbit
    void addOrbit(const Orbit &orbit)
    {
        _orbitList.push_back(orbit);
    }

    ///Returns number of orbits
    size_t size() const
    {
        return _orbitList.size();
    }

    /**
    Returns the number of orbits which are made up of N bodies
    */
    unsigned int getNumberOfNClusters(unsigned int N) const
    {
        unsigned int count = 0;
        for (const auto &orbit : _orbitList)
        {
            if (orbit.getRepresentativeCluster().getNumberOfBodies() == N)
            {
                count++;
            }
        }
        return count;
    }

    ///Return a copy of the orbit at position i in _orbitList
    Orbit getOrbit(unsigned int i) const
    {
        if (i >= size())
        {
            throw std::out_of_range("Error: Tried accessing orbit at out of bound index. Orbit OrbitList::getOrbit");
        }
        return _orbitList[i];
    }

    
    /// Clears the _orbitList
    void clear()
    {
        _orbitList.clear();
    }

    /// Sort the orbitlist
    void sort()
    {
        std::sort(_orbitList.begin(), _orbitList.end());
    }

    ///Returns orbitlist
    std::vector<Orbit> getOrbitList() const
    {
        return _orbitList;
    }

    int findOrbit(const Cluster &) const;

    /** 
    Prints information about the orbitlist


    */
    void print(int verbosity=0) const
    {
        int orbitCount =0;
        for(const auto &orbit : _orbitList)
        {
            std::cout<<"Orbit number: "<<orbitCount++<<std::endl;
            std::cout<<"Representative cluster "<<std::endl;
            orbit.getRepresentativeCluster().print();

            std::cout<<"Multiplicities: "<< orbit.size()<<std::endl;
            if(verbosity>1)
            {
                std::cout<<"Duplicates: "<< orbit.getNumberOfDuplicates()<<std::endl;
            }
            std::cout<<std::endl;
        }

    }
  private:
  int findOrbit(const Cluster &, const std::unordered_map<Cluster,int> &) const;

    std::vector<Orbit> _orbitList;
};