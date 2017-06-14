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

#include <boost/functional/hash.hpp>
using boost::hash;
using boost::hash_combine;
using boost::hash_value;

using namespace Eigen;

namespace py = pybind11;

class Cluster
{
  public:
    Cluster()
    {
        //empty constructor
    }

    Cluster(std::vector<int> &sites, std::vector<double> &distances)
    {
        _sites = sites;
        _distances = distances;
        sortCluster();
    }

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
        std::vector<std::pair<std::vector<std::pair<double, int>>, std::pair<int, int>>> first_dists(_sites.size());

        for (int i = 0; i < _sites.size(); i++)
        {
            first_dists[i] = std::make_pair(getDistsToSite(i), std::make_pair(_sites[i], i));
        }
        //ordering
        std::sort(first_dists.begin(), first_dists.end());

        std::vector<int> minimumOrder = getOrderFromFirstDists(first_dists[0]);
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

        //check if we have case 1
        //if not we do case 2 (both a check and a doer)

        auto equal_minimum_first_sites = getEqual_minimum_first_sites(first_dists);
        if (equal_minimum_first_sites.size() > 1)
        {
            std::cout << "case 1" << std::endl;
            auto case1_data = case1_min_indices(equal_minimum_first_sites);
            min_distance = std::get<1>(case1_data);
            min_sites = std::get<2>(case1_data);
        }
        else
        {
            std::cout << "case 2" << std::endl;
            ///Do case 2
            auto min_data = case2_min_indices(minimumOrder[0], first_dists[0].first);
            min_distance = std::get<1>(min_data);
            min_sites = std::get<2>(min_data);
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
        _distances = min_distance;
        _sites = min_sites;

        validateSorting();
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
    std::tuple<std::vector<int>, std::vector<double>, std::vector<int>> getMinimumStateBruteForce()
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

        auto minBruteForce = getMinimumStateBruteForce();
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
    std::vector<std::pair<std::vector<std::pair<double, int>>, std::pair<int, int>>> getEqual_minimum_first_sites(
        const std::vector<std::pair<std::vector<std::pair<double, int>>, std::pair<int, int>>> &first_dists) const
    {
        std::vector<std::pair<std::vector<std::pair<double, int>>, std::pair<int, int>>> equalFirstDists;
        equalFirstDists.push_back(first_dists[0]);
        int first_index_site = _sites[equalFirstDists[0].second.second];
        for (int i = 1; i < first_dists.size(); i++)
        {
            int this_site = _sites[first_dists[i].second.second];
            if (first_index_site == this_site && isEqualFirstDists(equalFirstDists[0].first, first_dists[i].first))
            {
                equalFirstDists.push_back(first_dists[i]);
            }
            else
            {
                break;
            }
        }

        std::cout<<" Found "<< equalFirstDists.size()<< " equal first dists"<<std::endl;
        std::cout<<"Minimum and second minimum: "<<std::endl;
        for(int i = 0; i < 2 ; i++)
        {   
            for(int j = 0; j <first_dists[i].first.size(); j++)
            {
                std::cout<<first_dists[i].first[j].first<< " "<<first_dists[i].first[j].second<< "  ";
            }
            std::cout<< _sites[first_dists[i].second.second]<<std::endl;
        }
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
        equal_minimum_first_dists
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

    std::tuple<std::vector<int>, std::vector<double>, std::vector<int>> case1_min_indices(const std::vector<std::pair<std::vector<std::pair<double, int>>, std::pair<int, int>>> &equal_minimum_first_dists) const
    {
        //Only one minimum dists..
        std::vector<std::pair<double, int>> minimum_dists;
        for (const auto &pair_index : equal_minimum_first_dists[0].first)
        {
            minimum_dists.push_back(std::make_pair(pair_index.first, _sites[pair_index.second]));
        }

        std::vector<double> min_distances;
        std::vector<int> min_sites;
        std::vector<int> min_order;
        if (isCase2(minimum_dists))
        {
            auto case_2_solution_first = case2_min_indices(equal_minimum_first_dists[0].second.second, equal_minimum_first_dists[0].first);

            min_order = std::get<0>(case_2_solution_first);
            min_distances = std::get<1>(case_2_solution_first);
            min_sites = std::get<2>(case_2_solution_first);
            //get first solution for comparing:
            //note that we begin on 1 since we did first manually
            for (int i = 1; i < equal_minimum_first_dists.size(); i++)
            {
                //case2_min_indices(const int i_index, const std::vector<std::pair<double, int>> &i_dist)
                auto case_2_solution = case2_min_indices(equal_minimum_first_dists[i].second.second, equal_minimum_first_dists[i].first);

                auto min_order_trial = std::get<0>(case_2_solution);
                auto min_distances_trial = std::get<1>(case_2_solution);
                auto min_sites_trial = std::get<2>(case_2_solution);
                if (compare_sites_dists(min_distances_trial, min_sites_trial, min_distances, min_sites))
                {
                    min_order = min_order_trial;
                    min_distances = min_distances_trial;
                    min_sites = min_sites_trial;
                }
            }
        }
        else
        {

            min_order = getOrderFromFirstDists(equal_minimum_first_dists[0]);
            min_distances = getReorderedDistances(min_order);
            min_sites = getReorderedSites(min_order);
            for (int i = 1; i < equal_minimum_first_dists.size(); i++)
            {
                auto min_order_trial = getOrderFromFirstDists(equal_minimum_first_dists[i]);
                {
                    auto min_order_trial = getOrderFromFirstDists(equal_minimum_first_dists[i]);
                    auto min_distances_trial = getReorderedDistances(min_order_trial);
                    auto min_sites_trial = getReorderedSites(min_order_trial);

                    if (compare_sites_dists(min_distances_trial, min_sites_trial, min_distances, min_sites))
                    {
                        min_order = min_order_trial;
                        min_distances = min_distances_trial;
                        min_sites = min_sites_trial;
                    }
                }
            }
        }
        return std::make_tuple(min_order, min_distances, min_sites);
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
    bool isCase2(const std::vector<std::pair<double, int>> &first_dists) const
    {
        return std::adjacent_find(first_dists.begin(), first_dists.end()) == first_dists.end();
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
    std::tuple<std::vector<int>, std::vector<double>, std::vector<int>> case2_min_indices(const int i_index, const std::vector<std::pair<double, int>> &i_dist) const
    {

        int counter;
        //this maps map distance, site to each indice in i_dist that have equal distances and sites
        std::map<std::pair<double, int>, std::vector<std::pair<int, int>>> uniqueDistsWithIndices;
        std::vector<int> minimumOrder;
        minimumOrder.reserve(_sites.size());
        minimumOrder.push_back(i_index); //this is always first for case2

        // create the current  mimumum order i_index , i_dist[0].second, ...
        for (const auto &dist_ind : i_dist)
        {
            minimumOrder.push_back(dist_ind.second);
        }
        std::vector<double> min_distances = getReorderedDistances(minimumOrder);
        std::vector<int> min_sites = getReorderedSites(minimumOrder);
        for (int i = 0; i < i_dist.size(); i++)
        {
            if (min_distances[i] != i_dist[i].first)
            {
                std::cout << "Error: min distances and i_dist do not equal" << std::endl;
                std::cout << min_distances[i] << " " << i_dist[i].first << std::endl;
                throw("Error: min distances and i_dist do not equal");
            }
        }

        //set up the uniqueDists
        // distance, site  = [(i_dist_index, i_nbr),... , ]
        counter = 0;
        for (const auto &dist_indice : i_dist)
        {
            uniqueDistsWithIndices[std::make_pair(dist_indice.first, _sites[dist_indice.second])].push_back(std::make_pair(counter++, dist_indice.second));
        }

        std::cout << "distances: " << std::endl;
        for (auto d : min_distances)
        {
            std::cout << d << " ";
        }
        std::cout << std::endl;
        std::cout << "i_nbrs: " << std::endl;
        for (auto d : i_dist)
        {
            std::cout << d.second << " ";
        }
        std::cout << std::endl;
        std::cout << "i_dist: " << std::endl;
        for (auto d : i_dist)
        {
            std::cout << d.first << " ";
        }
        std::cout << std::endl;

        //identical indices is a vector of vectors, first int is order it has in i_nbrs and second int is the index to the site in the cluster
        std::vector<std::vector<std::pair<int, int>>> identicalIndices;
        for (const auto &dist_indice_pair : uniqueDistsWithIndices)
        {
            std::cout << "( " << dist_indice_pair.first.first << " " << dist_indice_pair.first.second << " ) = [ ";
            for (auto id : dist_indice_pair.second)
            {
                std::cout << "(" << id.first << " " << id.second << " )";
            }
            std::cout << "] " << std::endl;
            if (dist_indice_pair.second.size() > 1)
            {
                identicalIndices.push_back(dist_indice_pair.second);
            }
        }
        //return if no dists, sites are equal
        if (identicalIndices.size() == 0)
        {
            // std::cout<< "no dists equal:" <<std::endl;
            return std::make_tuple(minimumOrder, min_distances, min_sites);
        }

        //sort identical indices so they appear in order in terms of indice in minimumOrder:
        //this is needed when taking the permutations of the identical sites
        for (auto &vec : identicalIndices)
        {
            std::sort(vec.begin(), vec.end());
        }
        //current minimum dists and sites

        //do all permutations of the identical dists, sites
        // identIndices = std::vector<std::pair<int, int>>
        for (auto &identIndices : identicalIndices)
        {

            //i_dist_indices are the sites that we can freely swap between
            std::vector<int> i_dist_indices;
            for (const auto &i_dist_index_i_nbr : identIndices)
            {
                // std::cout<<"( "<< i_dist_index_i_nbr.first<< " "<< i_dist_index_i_nbr.second<< ") ";
                i_dist_indices.push_back(i_dist_index_i_nbr.first);
            }
            // std::cout<<std::endl;
            do
            {
                std::vector<int> trial_indices = minimumOrder;
                for (int i = 0; i < i_dist_indices.size(); i++)
                {
                    trial_indices[i_dist_indices[i] + 1] = identIndices[i].second;
                }
                std::cout << "trial indices" << std::endl;
                for (auto i : trial_indices)
                {
                    std::cout << i << " ";
                }
                // std::cout<<":";
                // for(auto i : i_dist_indices){std::cout<<i<< "";}
                std::cout << std::endl;

                // for(const auto &i_dist_index_i_nbr : identIndices )
                // {
                //     trial_indices.push_back(i_dist_index_i_nbr.first);
                // }
                std::cout << "trial min dists " << std::endl;
                std::vector<double> trial_min_distances = getReorderedDistances(trial_indices);
                for (auto i : trial_min_distances)
                {
                    std::cout << i << " ";
                }
                std::cout << std::endl;
                std::vector<int> trial_min_sites = getReorderedSites(trial_indices);
                if (!compare_sites_dists(min_distances, min_sites, trial_min_distances, trial_min_sites))
                {
                    min_distances = trial_min_distances;
                    min_sites = trial_min_sites;
                }
            } while (std::next_permutation(i_dist_indices.begin(), i_dist_indices.end()));
        }

        return std::make_tuple(minimumOrder, min_distances, min_sites);
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
    Get all dists that origin from site i

    Returns a sorted double (distances) vector along with the index in _sites it points to
    */
    std::vector<std::pair<double, int>> getDistsToSite(int i)
    {
        std::vector<std::pair<double, int>> dists;
        int counter = 0;
        for (int k = 0; k < _sites.size(); k++)
        {
            for (int l = k + 1; l < _sites.size(); l++)
            {
                if (k == i or l == i)
                {
                    if (k != i)
                    {
                        dists.push_back(std::make_pair(_distances[counter], k));
                    }
                    else
                    {
                        dists.push_back(std::make_pair(_distances[counter], l));
                    }
                }
                counter++;
            }
        }
        if (counter != _distances.size())
        {
            std::cout << "Error: count not equal to distance size " << _distances.size() << " counter " << counter << std::endl;
            throw std::out_of_range("");
        }
        std::sort(dists.begin(), dists.end());
        return dists;
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

    // comparison operator for sortable clusters
    friend bool operator<(const Cluster &c1, const Cluster &c2)
    {
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
        std::cout << std::endl;
    }

  private:
    std::vector<int> _sites;
    std::vector<double> _distances;
    std::map<std::vector<int>, int> _element_counts;
    double symprec;
    double roundDouble(const double &double_value)
    {
        return round(double_value * 1.0 / symprec) / (1.0 / symprec);
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
