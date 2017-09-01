#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <vector>
#include "Orbit.hpp"
#include "ManybodyNeighborlist.hpp"
#include "Structure.hpp"
#include "Cluster.hpp"
#include "Neighborlist.hpp"
#include <unordered_map>
#include <unordered_set>
#include "LatticeNeighbor.hpp"
#include "hash_functions.hpp"
/**
Class OrbitList

contains a sorted vector or orbits


*/

class OrbitList
{
  public:
    OrbitList();
    OrbitList(const std::vector<Neighborlist> &neighborlists, const Structure &);
    OrbitList(const Structure &, const std::vector<std::vector<LatticeNeighbor>> &, const std::vector<Neighborlist> &);
    

    /**
    The strucute is a super cell
    The unsigned int is the index of the lattice site () (should probable be mapped to index 0 in  prim latnbr )
    the map maps primitive lattice neighbors to lattice neighbors in the supercell

    strategy is to get the translated orbit and then map it using the map and that should be the partial supercell orbit of this site
    add together all sites and you get the full supercell porbot
    */
    Orbit getSuperCellOrbit(const Structure &, const unsigned int, const std::unordered_map<LatticeNeighbor, LatticeNeighbor>&) const;

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
    void print(int verbosity = 0) const
    {
        int orbitCount = 0;
        for (const auto &orbit : _orbitList)
        {
            std::cout << "Orbit number: " << orbitCount++ << std::endl;
            std::cout << "Representative cluster " << std::endl;
            orbit.getRepresentativeCluster().print();

            std::cout << "Multiplicities: " << orbit.size() << std::endl;
            if (verbosity > 1)
            {
                std::cout << "Duplicates: " << orbit.getNumberOfDuplicates() << std::endl;
            }
            std::cout << std::endl;
        }
    }

    void addClusterToOrbitlist(const Cluster &cluster, const std::vector<LatticeNeighbor> &, std::unordered_map<Cluster, int> &);

    void addPermutationMatrixColumns(std::vector<std::vector<std::vector<LatticeNeighbor>>> &lattice_neighbors, std::unordered_set<std::vector<int>,VectorHash> &taken_rows, const std::vector<LatticeNeighbor> &lat_nbrs, const std::vector<int> &pm_rows,
                                     const std::vector<std::vector<LatticeNeighbor>> &permutation_matrix, const std::vector<LatticeNeighbor> &col1, bool) const;

    std::vector<LatticeNeighbor> getColumn1FromPM(const std::vector<std::vector<LatticeNeighbor>> &, bool sortIt = true) const;
    std::vector<int> findRowsFromCol1(const std::vector<LatticeNeighbor> &col1, const std::vector<LatticeNeighbor> &latNbrs, bool sortit = true) const;

    bool validatedCluster(const std::vector<LatticeNeighbor> &) const;
    void addOrbitsFromPM(const Structure &, const std::vector<std::vector<std::vector<LatticeNeighbor>>> &) ;
    void addOrbitFromPM(const Structure &,const std::vector<std::vector<LatticeNeighbor>> &);
    void checkEquivalentClusters(const Structure &) const;
    

    std::vector<LatticeNeighbor> translateSites(const std::vector<LatticeNeighbor> &, const unsigned int ) const;
    std::vector<std::vector<LatticeNeighbor>> getSitesTranslatedToUnitcell(const std::vector<LatticeNeighbor> &) const;
    std::vector<std::pair<std::vector<LatticeNeighbor>,std::vector<int>>> getMatchesInPM(const std::vector<std::vector<LatticeNeighbor>> &, const std::vector<LatticeNeighbor> &) const;

  private:
    int findOrbit(const Cluster &, const std::unordered_map<Cluster, int> &) const;
    Structure _primitiveStructure;
    std::vector<Orbit> _orbitList;
};