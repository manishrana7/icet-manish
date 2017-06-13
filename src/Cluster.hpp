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
        std::vector<int> minimumOrder(_sites.size());
        int min_index_count = 0;

        minimumOrder[min_index_count++] = first_dists[0].second.second;

        for (const auto &dist_pair : first_dists[0].first)
        {
            minimumOrder[min_index_count++] = dist_pair.second;
        }

        if (minimumOrder.size() != _sites.size())
        {
            throw("Error minimumorder.size != _sites.size()");
        }

        auto min_data = case2_min_indices(minimumOrder[0],first_dists[0].first);
        _distances = std::get<1>(min_data);
        _sites = std::get<2>(min_data);
        //setThisOrder(minimumOrder);
        if (_distances[0] != first_dists[0].first[0].first)
        {
            std::cout << "error distances not as expected, expected: "
                      << first_dists[0].first[0].first << ". Gotten: "
                      << _distances[0] << std::endl;
            throw("Error first dist is not the expected first dist");
        }
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
    std::tuple<std::vector<int>,std::vector<double>, std::vector<int>> case2_min_indices(const int i_index,const std::vector<std::pair<double, int>> &i_dist) const
    {

        int counter;
        //this maps map distance, site to each indice in i_dist that have equal distances and sites
        std::map<std::pair<double, int>, std::vector<std::pair<int, int>>> uniqueDistsWithIndices;
        std::vector<int> minimumOrder;
        minimumOrder.reserve(_sites.size());
        minimumOrder[0] = i_index; //this is always first for case2
        


        // create the current  mimumum order i_index , i_dist[0].second, ...
        for (const auto &dist_int : i_dist)
        {
            minimumOrder.push_back(dist_int.second);
        } 
        std::vector<double> min_distances = getReorderedDistances(minimumOrder);
        std::vector<int> min_sites = getReorderedSites(minimumOrder);

        //set up the uniqueDists
        // distance, site  = [(i_dist_index, i_nbr),... , ]
        counter = 0;
        for (const auto &dist_indice : i_dist)
        {
            uniqueDistsWithIndices[std::make_pair(dist_indice.first, _sites[dist_indice.second])].push_back(std::make_pair(counter++, dist_indice.second));        
        }

        std::cout<<"distances: "<<std::endl;
        for(auto d : min_distances){std::cout<<d<< " ";}
        std::cout<<std::endl;
        std::cout<<"i_nbrs: "<<std::endl;
        for(auto d : i_dist){std::cout<<d.second<< " ";}
        std::cout<<std::endl;
        std::cout<<"i_dist: "<<std::endl;
        for(auto d : i_dist){std::cout<<d.first<< " ";}
        std::cout<<std::endl;

        //identical indices is a vector of vectors, first int is order it has in i_nbrs and second int is the index to the site in the cluster        
        std::vector<std::vector<std::pair<int, int>>> identicalIndices;
        for (const auto &dist_indice_pair : uniqueDistsWithIndices)
        {
            std::cout<<"( "<<dist_indice_pair.first.first<< " "<<dist_indice_pair.first.second<<" ) = [ ";
            for(auto id : dist_indice_pair.second )
            {
                std::cout<<"("<<id.first<< " "<<id.second<< " )";
            }
            std::cout<<"] "<<std::endl;
            if (dist_indice_pair.second.size() > 1)
            {
                identicalIndices.push_back(dist_indice_pair.second);
            }
            
        }
        //return if no dists, sites are equal
        if (identicalIndices.size() == 0)
        {
            std::cout<< "no dists equal:" <<std::endl;
            return std::make_tuple(minimumOrder,min_distances, min_sites)  ;
        }

        //sort identical indices so they appear in order in terms of indice in minimumOrder:
        //this is needed when taking the permutations of the identical sites
        for(auto &vec : identicalIndices)
        {
            std::sort(vec.begin(), vec.end());
        }
        //current minimum dists and sites
        
        //do all permutations of the identical dists, sites
        // identIndices = std::vector<std::pair<int, int>>
        for(auto &identIndices : identicalIndices)
        {
            
            //i_dist_indices are the sites that we can freely swap between
            std::vector<int> i_dist_indices;
            for(const auto &i_dist_index_i_nbr : identIndices)
            {
                std::cout<<"( "<< i_dist_index_i_nbr.first<< " "<< i_dist_index_i_nbr.second<< ") ";
                i_dist_indices.push_back(i_dist_index_i_nbr.first);
            }
            std::cout<<std::endl;
            do {
                std::vector<int> trial_indices = minimumOrder;
                for(int i=0 ; i< i_dist_indices.size(); i++)
                {
                    trial_indices[i_dist_indices[i]+1] = identIndices[i].second;
                }
                for(auto i : trial_indices){std::cout<<i<< "";}
                std::cout<<":";
                for(auto i : i_dist_indices){std::cout<<i<< "";}
                std::cout<<std::endl;

                // for(const auto &i_dist_index_i_nbr : identIndices )
                // {
                //     trial_indices.push_back(i_dist_index_i_nbr.first);
                // }
                std::vector<double> trial_min_distances = getReorderedDistances(trial_indices);     
                std::vector<int> trial_min_sites = getReorderedSites(trial_indices);
                if( !compare_sites_dists(min_distances,min_sites, trial_min_distances, trial_min_sites ))
                {
                    min_distances = trial_min_distances;
                    min_sites = trial_min_sites;                    
                }
        } while ( std::next_permutation(i_dist_indices.begin(), i_dist_indices.end()) );

        }

        return std::make_tuple(minimumOrder,min_distances, min_sites);

    }


    /**

    compare distances and sites if dist1 and sites1 < dists2 and sites2
    (first compare distances, if they are equal compare sites)
    */

    bool compare_sites_dists(const std::vector<double> &dist1,const std::vector<int> &sites1,const std::vector<double> &dist2, const std::vector<int> &sites2 ) const
    {
        if( dist1 < dist2)
        {
            return true;
        }
        if( dist1 > dist2)
        {
            return false;
        }

        return sites1 < sites2;
    }

    /*
    Get all dists that origin from site i

    Returns a sorted double (distances) vector
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
