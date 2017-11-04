#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <pybind11/stl.h>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "Structure.hpp"
#include "LatticeNeighbor.hpp"
#include "Geometry.hpp"

#include <boost/functional/hash.hpp>
using boost::hash;
using boost::hash_combine;
using boost::hash_value;

using namespace Eigen;

namespace py = pybind11;

/**
Help struct which handle the I_Neighbors concept used in sorting a cluster.


*/
struct I_Neighbors
{

    I_Neighbors()
    {
        //empty constructor;
    }
    I_Neighbors(const int &this_index, const int &this_site, const std::vector<double> &distances, const std::vector<int> &sites, const std::vector<int> &indices)
    {
        i_index = this_index;
        i_site = this_site;

        _distances = distances;
        _sites = sites;
        _indices = indices;
    }

    std::vector<double> getIdists() const
    {
        return _distances;
    }

    std::vector<int> getSites() const
    {
        return _sites;
    }

    std::vector<int> getFullSites() const
    {
        std::vector<int> fullSites = {i_site};
        for (const auto &site : _sites)
        {
            fullSites.push_back(site);
        }
        return fullSites;
    }

    std::vector<int> getINeighbors() const
    {
        return _indices;
    }

    void setDistances(const std::vector<double> &distances)
    {
        _distances = distances;
    }

    void setSites(const std::vector<int> &sites)
    {
        _sites = sites;
    }
    void setIndices(const std::vector<int> &indices)
    {
        _indices = indices;
    }

    std::vector<int> getFullIndices() const
    {
        std::vector<int> fullIndices = {i_index};
        fullIndices.reserve(_indices.size());
        for (const auto &index : _indices)
        {
            fullIndices.push_back(index);
        }
        return fullIndices;
    }

    /**
    Looks for identical dists and sites

    returns a vector of vectors like:
    [i,j,k], [l, p] where i,j,k then have the same sites and distances
    and l,p have same sites and distances (but different from [i,j,k])

    The strategy is to use that the distances are sorted
    */
    std::vector<std::vector<int>> getIdenticalIndices() const
    {

        if (_distances.size() != _sites.size())
        {
            throw std::runtime_error("_distances and _sites not equal in size in getIdenticalIndices");
        }
        if (_indices.size() != _sites.size())
        {
            throw std::runtime_error("_distances and _sites not equal in size in getIdenticalIndices");
        }

        //maps dists, sites to vector of indices
        std::map<std::pair<double, int>, std::vector<int>> equalIndicesMap;
        for (int i = 0; i < _distances.size(); i++)
        {
            equalIndicesMap[std::make_pair(_distances[i], _sites[i])].push_back(i + 1);
        }

        std::vector<std::vector<int>> equalIndices;
        // if (equalIndicesMap.size() == _distances.size())
        // {
        //     //all dists, sites pairs were unique
        //     return equalIndices;
        // }

        for (const auto &mapPair : equalIndicesMap)
        {
            // std::cout << mapPair.first.first << " " << mapPair.first.second << " ";
            // for (auto d : mapPair.second)
            // {
            //     std::cout << d << " ";
            // }
            // std::cout << std::endl;
            if (mapPair.second.size() > 1)
            {
                equalIndices.push_back(mapPair.second);
            }
        }
        return equalIndices;
    }

    friend bool operator<(const I_Neighbors &i_n1, const I_Neighbors &i_n2)
    {
        if (i_n1._distances.size() < i_n2._distances.size())
        {
            return true;
        }
        else if (i_n1._distances.size() > i_n2._distances.size())
        {
            return false;
        }

        for (int i = 0; i < i_n1._distances.size(); i++)
        {

            // if (fabs(i_n1._distances[i] - i_n2._distances[i]) < 0.4)
            // {
            //     std::cout << i_n1._distances[i] << " < " << i_n2._distances[i] << " = " << std::boolalpha << (i_n1._distances[i] < i_n2._distances[i]) << std::endl;
            //     std::cout << i_n1._distances[i] << " > " << i_n2._distances[i] << " = " << std::boolalpha << (i_n1._distances[i] > i_n2._distances[i]) << std::endl;
            // }

            if (i_n1._distances[i] < i_n2._distances[i])
            {
                // std::cout << "return3" << std::endl;
                return true;
            }
            if (i_n1._distances[i] > i_n2._distances[i])
            {
                // std::cout << "return4" << std::endl;
                return false;
            }
        }

        if (i_n1.i_site < i_n2.i_site)
        {
            // std::cout << "return5" << std::endl;
            return true;
        }
        if (i_n1.i_site > i_n2.i_site)
        {
            // std::cout << "return6" << std::endl;
            return false;
        }

        for (int i = 0; i < i_n1._sites.size(); i++)
        {
            if (i_n1._sites[i] < i_n2._sites[i])
            {
                // std::cout << "return7" << std::endl;
                return true;
            }
            else if (i_n1._sites[i] > i_n2._sites[i])
            {
                // std::cout << "return8" << std::endl;
                return false;
            }
        }
        //everything is equal => return false;
        // std::cout << "return9" << std::endl;
        return false;
    }
    void print() const
    {

        for (auto d : _distances)
        {
            std::cout << d << " ";
        }
        std::cout << " : ";
        std::cout << i_site << " ";
        for (auto d : _sites)
        {
            std::cout << d << " ";
        }
        std::cout << i_index << " ";
        for (auto d : _indices)
        {
            std::cout << d << " ";
        }
        std::cout << std::endl;
    }

    friend bool operator==(const I_Neighbors &i_n1, const I_Neighbors &i_n2)
    {
        if (i_n1._distances != i_n2._distances)
        {
            return false;
        }
        else if (i_n1._sites != i_n2._sites)
        {
            return false;
        }
        return true;
    }

    std::vector<double> _distances;
    std::vector<int> _sites;
    std::vector<int> _indices;
    int i_index;
    int i_site;
};

class Cluster
{
  public:
    Cluster()
    {
        //empty constructor
    }

    // Cluster(std::vector<int> &sites, std::vector<double> &distances, const bool sortedCluster = true, const int clusterTag = 0)
    // {
    //     _symprec = 1e-5;
    //     _sites = sites;
    //     _distances = distances;
    //     _sortedCluster = sortedCluster;
    //     _clusterTag = clusterTag;
    //     if (_sortedCluster)
    //     {
    //         sortCluster();
    //     }

    //     //Run this if you doubt the sorting. It will bruteforce all
    //     //possible ways to see if there is another way to rearrange the cluster into a more "lower form"
    //     //validateSorting();
    // }

    ///Create cluster from a structure and latticeNeigbhors
    Cluster(const Structure &structure,
            const std::vector<LatticeNeighbor> &latticeNeighbors,
            const bool sortedCluster = true, const int clusterTag = 0);

    //counts the elements
    void count(const std::vector<int> &elements)
    {
        _element_counts[elements]++;
    }

    //get count of a specific element vector
    int getCount(const std::vector<int> &elements) const
    {
        const auto find = _element_counts.find(elements);
        if (find == _element_counts.end())
        {
            return 0;
        }
        else
        {
            return _element_counts.at(elements);
        }
    }
    //return the unique sites
    std::vector<int> getSites() const
    {
        return _sites;
    }

    //return the cluster distances
    std::vector<double> getDistances() const
    {
        return _distances;
    }

    /**
    Sorts the clusters internal sites and corresponding distances 
    so it is in the lowest form possible.

    This enables sorting and uniqueness of clusters since you can easily compare it using
    the defined compare and equal function.

    Algorithm for sorting:
    ----------------------
    First each site in the cluster gets its distances in sorted order to the other sites,
    denote this by i_dists, i_nbrs.

    if there is one smallest i_dist then i is the first site and i_nbrs[0] is the 
    second site and so on

    There are then some special cases which need to be adhered:

    case 1: 
        Two or more i_dists are equal

        Problem: 
            Which site is then first?

    case 2:
        Two or more distances are equal in i_dists.

        Problem:
            unclear which site is site 2,3 (if i_dist[0] == i_dist[1] )

    case 3:
        two or more distances are equal in j_dist where j is in i_nbrs

        Problem:
        What is the problem here?

    solutions case 1:
        1) assume case2 and case3 is not active.
            a) get "getDistIndices" and getReorderedSites for each case
            b) compare to each other. The lowest should be the one to use
        2) if case2 is active also:
           a )for each smallest i_dist:
                get solution from case 2
            b) compare case2 solutions of i_dist to find minimum.

    solutions case 2:
        for each j,k,l.. indices in i_nbr that have i_dists equal:
            for combination in all_combinations(j,k,l,...):
                get "getDistIndices" and getReorderedSites
            take the get_dist_indices and reordered sites that are smallest and 
            use the corresponding combination of indices.


    */
    void sortCluster()
    {
        if (_distances.size() == 0)
        {
            return;
        }
        std::vector<I_Neighbors> first_dists(_sites.size());

        for (int i = 0; i < _sites.size(); i++)
        {
            first_dists[i] = getDistsToSite(i);
        }
        //ordering
        // std::cout << "before sort" << std::endl;
        // for (auto d : first_dists)
        // {
        //     d.print();
        // }

        std::sort(first_dists.begin(), first_dists.end());

        // std::cout << "after sort" << std::endl;
        // for (auto d : first_dists)
        // {
        //     d.print();
        // }

        // first_dists[0].print();
        // first_dists[1].print();
        // std::cout << "ok test these two above " << std::endl;
        // std::cout << " less than? " << std::boolalpha << (first_dists[0] < first_dists[1]) << std::endl;
        std::vector<int> minimumOrder = first_dists[0].getFullIndices(); // getOrderFromFirstDists(first_dists[0]);
        // int min_index_count = 0;

        // minimumOrder[min_index_count++] = first_dists[0].second.second;

        // for (const auto &dist_pair : first_dists[0].first)
        // {
        //     minimumOrder[min_index_count++] = dist_pair.second;
        // }

        if (minimumOrder.size() != _sites.size())
        {
            throw("Error minimumorder.size != _sites.size()");
        }

        std::vector<double> min_distance;
        std::vector<int> min_sites;
        std::vector<int> min_indices;

        //check if we have case 1
        //if not we do case 2 (both a check and a doer)

        auto equal_minimum_first_sites = getEqual_minimum_first_sites(first_dists);
        if (equal_minimum_first_sites.size() > 1)
        {
            // std::cout << "case 1" << std::endl;
            auto min_data = case1_min_indices(equal_minimum_first_sites);
            min_distance = std::get<0>(min_data);
            min_sites = std::get<1>(min_data);
            min_indices = std::get<2>(min_data);
        }
        else
        {
            // std::cout << "case 2" << std::endl;
            ///Do case 2
            auto min_data = case2_min_indices(first_dists[0]);
            min_distance = std::get<0>(min_data);
            min_sites = std::get<1>(min_data);
            min_indices = std::get<2>(min_data);
        }

        /// some validation of algorithm and debugging
        // auto original_distance_copy = _distances;
        // std::sort(min_distance.begin(), min_distance.end());
        // std::sort(original_distance_copy.begin(), original_distance_copy.end());
        // for (int i = 0; i < original_distance_copy.size(); i++)
        // {
        //     if (original_distance_copy[i] != min_distance[i])
        //     {
        //         std::cout << "error: the distances havee been changed (Not only order)" << std::endl;
        //         throw std::runtime_error("error: the distances havee been changed (Not only order)");
        //     }
        // }
        if (min_distance.size() != _distances.size())
        {
            std::cout << "min_distances.size() = " << min_distance.size() << " _distances.size() = " << _distances.size() << std::endl;
            throw std::runtime_error("Error: algorithm ended up with wrong distance size");
        }

        if (min_sites.size() != _sites.size())
        {
            throw std::runtime_error("Error: algorithm ended up with wrong sites size");
        }
        _distances = min_distance;
        _sites = min_sites;
        // std::cout << " sort cluster found order: " << std::endl;
        // for (auto d : min_indices)
        // {
        //     std::cout << d << " ";
        // }
        // std::cout << std::endl;
        //setThisOrder(minimumOrder);
        // for (int i = 0; i < first_dists[0].first.size(); i++)
        // {
        //     if (_distances[i] != first_dists[0].first[i].first)
        //     {
        //         std::cout << "error distances not as expected, expected: "
        //                   << first_dists[0].first[i].first << ". Gotten: "
        //                   << _distances[i] << std::endl;
        //         throw std::runtime_error("Error first dist is not the expected first dist");
        //     }
        // }
    }

    /**
    Bruteforce attempt to find minumum cluster
    */
    std::tuple<std::vector<int>, std::vector<double>, std::vector<int>> getNumberOfAllowedComponentsnimumStateBruteForce()
    {

        std::vector<int> atomic_order(_sites.size());
        for (int i = 0; i < atomic_order.size(); i++)
        {
            atomic_order[i] = i;
        }
        auto min_order = atomic_order;
        auto min_dists = getReorderedDistances(atomic_order);
        auto min_sites = getReorderedSites(atomic_order);

        do
        {
            auto dists = getReorderedDistances(atomic_order);
            auto sites = getReorderedSites(atomic_order);
            if (compare_sites_dists(dists, sites, min_dists, min_sites))
            {
                min_dists = dists;
                min_sites = sites;
                min_order = atomic_order;
            }

        } while (std::next_permutation(atomic_order.begin(), atomic_order.end()));

        return std::make_tuple(min_order, min_dists, min_sites);
    }

    /**
    This validates the sorting by testing all possible 
    combinations of sites and sees if there is a lower state of the cluster.

    */
    void validateSorting()
    {

        auto minBruteForce = getNumberOfAllowedComponentsnimumStateBruteForce();
        auto bruteforce_order = std::get<0>(minBruteForce);
        auto bruteForceDists = std::get<1>(minBruteForce);
        auto bruteForceSites = std::get<2>(minBruteForce);

        if (compare_sites_dists(bruteForceDists, bruteForceSites, _distances, _sites))
        {
            std::cout << " bruteforce dists, sites and order:" << std::endl;
            for (auto d : bruteForceDists)
            {
                std::cout << d << " ";
            }
            std::cout << std::endl;
            for (auto d : bruteForceSites)
            {
                std::cout << d << " ";
            }
            std::cout << std::endl;
            for (auto d : bruteforce_order)
            {
                std::cout << d << " ";
            }
            std::cout << std::endl;
            std::cout << " algorithm result for dists and sites:" << std::endl;
            for (auto d : _distances)
            {
                std::cout << d << " ";
            }
            std::cout << std::endl;
            for (auto d : _sites)
            {
                std::cout << d << " ";
            }
            std::cout << std::endl;

            throw std::runtime_error("Error: brute force found a smaller");
        }
    }

    /**
    Get all equal first_dists that are identical to the minimum one
    By equality one mean that all distances and sites are equalr
    
    */
    std::vector<I_Neighbors> getEqual_minimum_first_sites(
        const std::vector<I_Neighbors> &i_neighbors) const
    {
        std::vector<I_Neighbors> equalFirstDists;
        equalFirstDists.push_back(i_neighbors[0]);

        for (int i = 1; i < i_neighbors.size(); i++)
        {

            if (i_neighbors[i] == i_neighbors[0])
            {
                equalFirstDists.push_back(i_neighbors[i]);
            }
            else
            {
                break;
            }
        }

        // std::cout << " Found " << equalFirstDists.size() << " equal first dists" << std::endl;
        // std::cout << "Minimum and second minimum: " << std::endl;
        // for (int i = 0; i < i_neighbors.size(); i++)
        // {
        //     i_neighbors[i].print();
        // }

        return equalFirstDists;
    }

    /**
   Returns true if dist_index1 is equal to dist_index2

   THe hitch is that dist_index1[i].second is a index in the cluster and has to be checked what site is on that
   index when comparing

   */
    bool isEqualFirstDists(const std::vector<std::pair<double, int>> &dist_index1, const std::vector<std::pair<double, int>> &dist_index2) const
    {
        for (int i = 0; i < dist_index1.size(); i++)
        {
            if (!(dist_index1[i].first == dist_index2[i].first && _sites[dist_index1[i].second] == _sites[dist_index2[i].second]))
            {
                return false;
            }
        }
        return true;
    }

    /**


    Arguments:
        equal_minimum_i_neighbors
        a vector of all the minimum first dists
    Do case 1: 
        solutions case 1:
        1) assume case2 and case3 is not active.
            a) get "getDistIndices" and getReorderedSites for each case
            b) compare to each other. The lowest should be the one to use
        2) if case2 is active also:
           a )for each smallest i_dist:
                get solution from case 2
            b) compare case2 solutions of i_dist to find minimum.

    */

    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> case1_min_indices(const std::vector<I_Neighbors> &equal_minimum_i_neighbors) const
    {

        auto minOrder = equal_minimum_i_neighbors[0].getFullIndices();
        auto minDistance = getReorderedDistances(minOrder);
        auto minSites = getReorderedSites(minOrder);

        if (isCase2(equal_minimum_i_neighbors[0]))
        {
            // std::cout << "case 2 inside case1" << std::endl;
            for (int i = 0; i < equal_minimum_i_neighbors.size(); i++)
            {
                //case2_min_indices(const int i_index, const std::vector<std::pair<double, int>> &i_dist)
                auto case_2_solution = case2_min_indices(equal_minimum_i_neighbors[i]);

                auto min_distances_trial = std::get<0>(case_2_solution);
                auto min_sites_trial = std::get<1>(case_2_solution);
                auto min_order_trial = std::get<2>(case_2_solution);

                if (compare_sites_dists(min_distances_trial, min_sites_trial, minDistance, minSites))
                {
                    minOrder = min_order_trial;
                    minDistance = min_distances_trial;
                    minSites = min_sites_trial;
                }
            }
        }
        else
        {
            // std::cout << "NOT case 2 inside case1" << std::endl;
            for (int i = 0; i < equal_minimum_i_neighbors.size(); i++)
            {

                auto min_order_trial = equal_minimum_i_neighbors[i].getFullIndices();
                auto min_distances_trial = getReorderedDistances(min_order_trial);
                auto min_sites_trial = getReorderedSites(min_order_trial);

                if (compare_sites_dists(min_distances_trial, min_sites_trial, minDistance, minSites))
                {
                    minOrder = min_order_trial;
                    minDistance = min_distances_trial;
                    minSites = min_sites_trial;
                }
            }
        }
        return std::make_tuple(minDistance, minSites, minOrder);
    }

    /**
    Gets the min order of indices as given by this "first dists"

    this format is used throughout this code to simpler extract the indices from a given first dist (i_dists)

    The format is 
    first_dists.first = a vector of distances and indices to the corresponding site the distance is to
    first_dists.second.first = the site of this index
    first_dists.second.second the index of this first_dists, i.e. the i in i_dists
    */

    std::vector<int> getOrderFromFirstDists(const std::pair<std::vector<std::pair<double, int>>, std::pair<int, int>> &first_dists) const
    {
        std::vector<int> minOrder(_sites.size());

        int counter = 0;
        minOrder[counter++] = first_dists.second.second;
        for (const auto &dist_index_pair : first_dists.first)
        {
            minOrder[counter++] = dist_index_pair.second;
        }
        return minOrder;
    }
    /** 
        Checks if these distance, sites vector has some equal elements 
         and is thus subject to case2
    */
    bool isCase2(const I_Neighbors &i_neighbors) const
    {
        for (int i = 0; i < i_neighbors._distances.size() - 1; i++)
        {
            if (i_neighbors._distances[i] == i_neighbors._distances[i + 1] && i_neighbors._sites[i] == i_neighbors._sites[i + 1])
            {
                return true;
            }
        }
        return false;
        //return std::adjacent_find(first_dists.begin(), first_dists.end()) == first_dists.end();
    }
    /**
    Rearranges distances and sites to the order given in minimumOrder

    example:
    -------
    minimumOrder[0] = current index which should be first place
    minimumOrder[1] = current index which should be second place
    and so on ... 

    This is a slightly non trivial operation since swapping two sites changes 
    the order of the distances in a complex way.

    The distances are ordered by:
     d_12, d_13, d_14, ... , d_23, d_24, ...


     solution:
     ---------
     The algorithm creates a vector of tuples with i, j, distance

     the current state is: 1, 2, d_12
                           1, 3, d_13

    
    which is then mapped to:
                           minOrder.index_of(1), minOrder.index_of(2), d_12
                           minOrder.index_of(1), minOrder.index_of(3), d_13


    where minOrder.index_of(1) = index of 1 in minOrder.
    this vector is then sorted and the positions will align to the new order.

    */
    void setThisOrder(const std::vector<int> &minimumOrder)
    {
        _sites = getReorderedSites(minimumOrder);
        _distances = getReorderedDistances(minimumOrder);
    }

    /**
    The distances are ordered by:

     d_12, d_13, d_14, ... , d_23, d_24, ...
     
     return a tuple of the current state:
                           1, 2, d_12
                           1, 3, d_13

    */
    std::vector<std::tuple<int, int, double>> getDistIndices() const
    {

        std::vector<std::tuple<int, int, double>> dist_indices;
        int counter = 0;
        for (int k = 0; k < _sites.size(); k++)
        {
            for (int l = k + 1; l < _sites.size(); l++)
            {
                dist_indices.push_back(std::make_tuple(k, l, _distances[counter++]));
            }
        }
        return dist_indices;
    }

    /**
     The algorithm creates a vector of tuples with i, j, distance

     the current state is: 1, 2, d_12
                           1, 3, d_13

    
    which is then mapped to:
                           indiceOrder.index_of(1), indiceOrder.index_of(2), d_12
                           indiceOrder.index_of(1), indiceOrder.index_of(3), d_13


    where indiceOrder.index_of(1) = index of 1 in indiceOrder.
    this vector is then sorted and the positions will align to the new order.


    */
    std::vector<std::tuple<int, int, double>> getDistIndices(const std::vector<int> &indiceOrder) const
    {

        std::vector<std::tuple<int, int, double>> dist_indices;
        dist_indices.reserve(_distances.size());
        int counter = 0;
        for (int k = 0; k < _sites.size(); k++)
        {
            for (int l = k + 1; l < _sites.size(); l++)
            {

                auto find_first = std::find(indiceOrder.begin(), indiceOrder.end(), k);
                auto find_second = std::find(indiceOrder.begin(), indiceOrder.end(), l);

                int first = std::distance(indiceOrder.begin(), find_first);
                int second = std::distance(indiceOrder.begin(), find_second);

                if (second < first)
                {
                    auto temp = first;
                    first = second;
                    second = temp;
                }
                dist_indices.push_back(std::make_tuple(first, second, _distances[counter++]));
            }
        }

        std::sort(dist_indices.begin(), dist_indices.end());

        return dist_indices;
    }

    /**
    Get the distances if the sites would have been rearranged according to indiceOrder

    returns a double (distance ) vector
    */

    std::vector<double> getReorderedDistances(const std::vector<int> &indiceOrder) const
    {
        auto dist_indices = getDistIndices(indiceOrder);
        int counter = 0;
        std::vector<double> reorderedDistances(_distances.size());
        for (const auto &tup : dist_indices)
        {
            reorderedDistances[counter++] = std::get<2>(tup);
        }
        return reorderedDistances;
    }

    /**
    Return the sites if the order would have been as is given in indiceOrder

    */
    std::vector<int> getReorderedSites(const std::vector<int> &indiceOrder) const
    {
        std::vector<int> tempSites(_sites.size());
        for (int i = 0; i < indiceOrder.size(); i++)
        {
            tempSites[i] = _sites[indiceOrder[i]];
        }
        return tempSites;
    }

    /**
    solutions case 2:
        for each j,k,l.. indices in i_nbr that have i_dists equal:
            for combination in all_combinations(j,k,l,...):
                get "getDistIndices" and getReorderedSites
            take the get_dist_indices and reordered sites that are smallest and 
            use the corresponding combination of indices.

    returns min_indices, correspong distances and corresponding sites
    */
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> case2_min_indices(const I_Neighbors &i_neighbor) const
    {

        std::vector<int> minimumOrder = i_neighbor.getFullIndices();

        std::vector<double> min_distances = getReorderedDistances(minimumOrder);
        std::vector<int> min_sites = getReorderedSites(minimumOrder);
        std::vector<int> min_indices = minimumOrder;
        // std::cout << "original distances " << std::endl;
        // for (auto d : min_distances)
        // {
        //     std::cout << d << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "printing i_neighbor " << std::endl;
        // i_neighbor.print();
        // std::cout << " Min order" << std::endl;
        // for (auto d : minimumOrder)
        // {
        //     std::cout << d << " ";
        // }
        // std::cout << std::endl;

        //identical indices is a vector of vectors
        // the identical indices are given in the global indices in the cluster
        std::vector<std::vector<int>> identicalIndices = i_neighbor.getIdenticalIndices();

        //return if no dists, sites are equal
        if (identicalIndices.size() == 0)
        {
            // std::cout << "no dists equal:" << std::endl;
            // i_neighbor.print();
            return std::make_tuple(min_distances, min_sites, min_indices);
        }

        // for (auto identVec : identicalIndices)
        // {
        //     for (auto d : identVec)
        //     {
        //         std::cout << d << " ";
        //     }
        //     std::cout << std::endl;
        // }

        if (minimumOrder.size() != _sites.size())
        {
            throw std::runtime_error("Minimumorder in case 2 is not equal in size to _sites");
        }
        //sort identical indices so they appear in order in terms of indice in minimumOrder:
        //this is needed when taking the permutations of the identical sites
        for (auto &vec : identicalIndices)
        {
            std::sort(vec.begin(), vec.end());
        }

        //do all permutations of the identical dists, sites
        // identIndices = std::vector<std::pair<int, int>>

        //this i_neighbor has the right first indice and site

        auto trial_order = minimumOrder;
        findMinimumIndicePermutation(0, identicalIndices, minimumOrder, trial_order, min_distances, min_sites, min_indices);

        return std::make_tuple(min_distances, min_sites, min_indices);
    }

    void findMinimumIndicePermutation(int currentIndiceSet,
                                      const std::vector<std::vector<int>> &identicalIndices,
                                      const std::vector<int> &minimumOrder,
                                      std::vector<int> trial_order,
                                      std::vector<double> &min_distances,
                                      std::vector<int> &min_sites,
                                      std::vector<int> &min_indices) const
    {
        std::vector<int> identicalIndiceSet = identicalIndices[currentIndiceSet];
        //trial_order = minimumOrder;
        do
        {
            // trial_order = minimumOrder;
            for (int i = 0; i < identicalIndiceSet.size(); i++)
            {
                trial_order[identicalIndiceSet[i]] = minimumOrder[identicalIndices[currentIndiceSet][i]];
            }

            auto distances_trial = getReorderedDistances(trial_order);
            auto trial_sites = getReorderedSites(trial_order);

            if (compare_sites_dists(distances_trial, trial_sites, min_distances, min_sites))
            {
                min_distances = distances_trial;
                min_sites = trial_sites;
                min_indices = trial_order;
            }
            if (currentIndiceSet + 1 < identicalIndices.size())
            {
                findMinimumIndicePermutation(currentIndiceSet + 1, identicalIndices, minimumOrder, trial_order, min_distances, min_sites, min_indices);
            }

        } while (std::next_permutation(identicalIndiceSet.begin(), identicalIndiceSet.end()));
    }
    /**

    compare distances and sites if dist1 and sites1 < dists2 and sites2
    (first compare distances, if they are equal compare sites)
    */

    bool compare_sites_dists(const std::vector<double> &dist1, const std::vector<int> &sites1, const std::vector<double> &dist2, const std::vector<int> &sites2) const
    {

        for (int i = 0; i < dist1.size(); i++)
        {
            if (dist1[i] < dist2[i])
            {
                // std::cout<<"true "<< dist1[i]<< " "<< dist2[i]<<std::endl;
                return true;
            }
            else if (dist1[i] > dist2[i])
            {
                // std::cout<<"false "<< dist1[i]<< " "<< dist2[i]<<std::endl;
                return false;
            }
        }
        if (dist1 < dist2)
        {
            // std::cout<<"true: dist1< dist2"<<std::endl;
            return true;
        }
        if (dist1 > dist2)
        {
            // std::cout<<"false: dist1 > dist2"<<std::endl;
            return false;
        }

        // std::cout<<"return sites sites1 < sites2 "<< std::boolalpha<< (sites1 < sites2)<<std::endl;

        return sites1 < sites2;
    }

    /*
    Get all dists, sites and indices that origin from site i

    Returns a I_Neigbhors structure
    */
    I_Neighbors getDistsToSite(const int i_index)
    {
        std::vector<std::tuple<double, int, int>> dists_site_index;
        int counter = 0;
        for (int k = 0; k < _sites.size(); k++)
        {
            for (int l = k + 1; l < _sites.size(); l++)
            {
                if (k == i_index or l == i_index)
                {
                    if (k != i_index)
                    {
                        dists_site_index.push_back(std::make_tuple(_distances[counter], _sites[k], k)); //Change to sites here or else it will mess up the sorting
                    }
                    else
                    {
                        dists_site_index.push_back(std::make_tuple(_distances[counter], _sites[l], l));
                    }
                }
                counter++;
            }
        }
        if (counter != _distances.size())
        {
            std::cout << "Error: count not equal to distance size " << _distances.size() << " counter " << counter << std::endl;
            std::cout << "Sites size: " << _sites.size() << std::endl;
            throw std::out_of_range(" count not equal to distance size");
        }
        std::sort(dists_site_index.begin(), dists_site_index.end());
        std::vector<double> distances(dists_site_index.size());
        std::vector<int> sites(dists_site_index.size());
        std::vector<int> indices(dists_site_index.size());

        for (int i = 0; i < dists_site_index.size(); i++)
        {
            distances[i] = std::get<0>(dists_site_index[i]);
            sites[i] = std::get<1>(dists_site_index[i]);
            indices[i] = std::get<2>(dists_site_index[i]);
        }
        I_Neighbors i_nbrs = I_Neighbors(i_index, _sites[i_index], distances, sites, indices);

        return i_nbrs;
    }

    /*
    Swap sites of site i and site j and also changes the corresponding distances.
    */
    void swapSites(const int i, const int j)
    {
        //swap sites
        std::iter_swap(_sites.begin() + i, _sites.begin() + j);

        std::vector<std::tuple<int, int, double>> dist_indices;
        int counter = 0;
        for (int k; k < _sites.size(); k++)
        {
            for (int l = k; l < _sites.size(); l++)
            {
                int first = k;
                int second = l;
                if (first == i)
                {
                    first = j;
                }
                if (first == j)
                {
                    first = i;
                }
                if (second == i)
                {
                    second = j;
                }
                if (second == j)
                {
                    second = i;
                }

                if (second < first)
                {
                    auto temp = first;
                    first = second;
                    second = temp;
                }
                dist_indices.push_back(std::make_tuple(first, second, _distances[counter++]));
            }
        }

        std::sort(dist_indices.begin(), dist_indices.end());
        counter = 0;
        for (const auto &tup : dist_indices)
        {
            _distances[counter++] = std::get<2>(tup);
        }
    }

    friend bool operator==(const Cluster &c1, const Cluster &c2)
    {
        if (c1.isSorted() != c2.isSorted())
        {
            throw std::runtime_error("Undefined behavior: comparing sorted and non-sorted cluster");
        }

        //Non-sorted clusters uses clustertag for comparison
        if (!c1.isSorted())
        {
            return c1.getClusterTag() == c2.getClusterTag();
        }

        ///@TODO: Think if this is enough
        if (c1._sites != c2._sites)
        {
            return false;
        }
        if (c1._distances != c2._distances)
        {
            return false;
        }
        return true;

        //return( c1 < c2 && c2 < c1);
    }

    friend bool operator!=(const Cluster &c1, const Cluster &c2)
    {
        return !(c1 == c2);
    }

    // comparison operator for sortable clusters
    friend bool operator<(const Cluster &c1, const Cluster &c2)
    {

        if (c1.isSorted() != c2.isSorted())
        {
            throw std::runtime_error("Undefined behavior: comparing sorted and non-sorted cluster");
        }

        if (!c1.isSorted())
        {
            return c1.getClusterTag() < c2.getClusterTag();
        }

        //1) compare number of bodies in cluster
        if (c1._sites.size() < c2._sites.size())
        {
            return true;
        }
        if (c1._sites.size() > c2._sites.size())
        {
            return false;
        }

        //2) compare distances
        for (int i = 0; i < c1._distances.size(); i++)
        {
            if (c1._distances[i] < c2._distances[i])
            {
                return true;
            }
            if (c1._distances[i] > c2._distances[i])
            {
                return false;
            }
        }

        //3) compare sites
        for (int i = 0; i < c1._sites.size(); i++)
        {
            if (c1._sites[i] < c2._sites[i])
            {
                return true;
            }
            if (c1._sites[i] > c2._sites[i])
            {
                return false;
            }
        }

        //4) if we are here then everything is equal so return false
        return false;
    }

    /**
    Print the cluster to standard out.

    Format is first distances then sites
    */
    void print() const
    {
        for (const auto d : _distances)
        {
            std::cout << d << " ";
        }
        std::cout << " :: ";
        for (const auto s : _sites)
        {
            std::cout << s << " ";
        }
        std::cout<< _geometricalSize<< " ";
        std::cout << std::endl;

    }

    ///Return true if this is a sorted cluster
    bool isSorted() const
    {
        return _sortedCluster;
    }

    ///Return the cluster tag used for identification if not dists/sites is used to distinguish
    int getClusterTag() const
    {
        return _clusterTag;
    }
    ///Sets the clustertag and marks this cluster as non-sorted
    void setClusterTag(const int clusterTag) 
    {
        _sortedCluster = false;
        _clusterTag = clusterTag;
    }
    ///Return the number of bodies (size of sites) of the cluster
    unsigned int getNumberOfBodies() const
    {
        return _sites.size();
    }

    double getGeometricalSize() const
    {
        return _geometricalSize;
    }

  private:
    std::vector<int> _sites;
    std::vector<double> _distances;
    std::map<std::vector<int>, int> _element_counts;
    double _geometricalSize;
    bool _sortedCluster;
    int _clusterTag;
    double _symprec;
    double roundDouble(const double &double_value)
    {
        return round(double_value * 1.0 / _symprec) / (1.0 / _symprec);
    }
};

namespace std
{
template <>
struct hash<Cluster>
{
    size_t
    operator()(const Cluster &k) const
    {

        // Compute individual hash values for first,
        // second and third and combine them using XOR
        // and bit shifting:
        size_t seed = 0;

        //if unsorted just use the clusterTag as seed
        if (!k.isSorted())
        {
            hash_combine(seed, k.getClusterTag());
            return seed;
        }

        for (const auto &distance : k.getDistances())
        {
            hash_combine(seed, hash_value(distance));
        }
        for (const auto &site : k.getSites())
        {
            hash_combine(seed, hash_value(site));
        }
        return seed;
    }
};
}
