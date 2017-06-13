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
    
    */
    void sortCluster()
    {
        std::vector<std::pair<std::vector<std::pair<double,int>>, std::pair<int, int>>> first_dists;

        for (int i = 0; i < _sites.size(); i++)
        {
            first_dists.push_back(std::make_pair(getDistsToSite(i), std::make_pair(_sites[i], i)));
        }
        //ordering
        std::sort(first_dists.begin(), first_dists.end());
        std::vector<int> minimumOrder(_sites.size());
        int min_index_count = 0;


        minimumOrder[min_index_count++] = first_dists[0].second.second;
        
        for(const auto &dist_pair : first_dists[0].first)
        {
            minimumOrder[min_index_count++] = dist_pair.second;        
        }
        

        if(minimumOrder.size() != _sites.size())
        {
            throw("Error minimumorder.size != _sites.size()");
        }

        setThisOrder(minimumOrder);
        if( _distances[0] != first_dists[0].first[0].first)
        {
             std::cout<<"error distances not as expected, expected: "
                      <<first_dists[0].first[0].first<< ". Gotten: "
                      << _distances[0]<<std::endl;
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
        //swap sites        
        std::vector<int> tempSites(_sites.size());
        for(int i = 0; i < minimumOrder.size(); i++)
        {
            tempSites[i] = _sites[minimumOrder[i]];
        }
        _sites = tempSites;


        std::vector<std::tuple<int, int, double>> dist_indices;
        int counter = 0;
        for (int k = 0; k < _sites.size(); k++)
        {
            for (int l = k + 1; l < _sites.size(); l++)
            {

                auto find_first  = std::find(minimumOrder.begin(), minimumOrder.end(), k);
                auto find_second = std::find(minimumOrder.begin(), minimumOrder.end(), l);

                int first =  std::distance(minimumOrder.begin(),  find_first);
                int second = std::distance(minimumOrder.begin(), find_second);

                if( second < first )
                {
                    auto temp = first;
                    first = second;
                    second = temp;
                }
                dist_indices.push_back(std::make_tuple(first, second, _distances[counter++]));
            }
        }
        //debugging cout
        // std::cout<<"set this order before: "<<std::endl;
        // for(const auto tup : dist_indices)
        // {
        //     std::cout<<std::get<0>(tup)<< " "<<std::get<1>(tup)<< " "<<std::get<2>(tup)<<std::endl;
        // }


        //sort to get the new order of dists
        std::sort(dist_indices.begin(), dist_indices.end());

        //debugging cout
        // std::cout<<"set this order after sort: "<<std::endl;
        // for(const auto tup : dist_indices)
        // {
        //     std::cout<<std::get<0>(tup)<< " "<<std::get<1>(tup)<< " "<<std::get<2>(tup)<<std::endl;
        // }

        counter = 0;
        for (const auto &tup : dist_indices)
        {
            _distances[counter++] = std::get<2>(tup);
        }
    }

    /*
    Get all dists that origin from site i

    Returns a sorted double (distances) vector
    */
    std::vector<std::pair<double,int>> getDistsToSite(int i)
    {
        std::vector<std::pair<double,int>> dists;
        int counter = 0;
        for (int k=0; k < _sites.size(); k++)
        {
            for (int l = k+1; l < _sites.size(); l++)
            {
                if (k == i or l == i)
                {
                    if(k != i)
                    {
                        dists.push_back(std::make_pair(_distances[counter],k));
                    }
                    else
                    {
                        dists.push_back(std::make_pair(_distances[counter],l));
                    }
                    
                }
                counter++;
            }
        }
        if ( counter != _distances.size())
        {
            std::cout<<"Error: count not equal to distance size "<< _distances.size() << " counter "<<counter<<std::endl;
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
