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

    void count(const std::vector<int> &elements)
    {
        _element_counts[elements]++;
    }

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

    std::vector<int> getSites() const
    {
        return _sites;
    }

    std::vector<double> getDistances() const
    {
        return _distances;
    }

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

        //minimumOrder[min_index_count++] = first_dists[0].first[0].second; //first min index
        // first_dists[0] = std::pair<std::vector<std::pair<double,int>>, std::pair<int, int>>
        // first_dists[0].first = std::vector<std::pair<double,int>

        minimumOrder[min_index_count++] = first_dists[0].second.second;
        // std::cout<<" distance,int ("<<first_dists[0].second.second<< ")";
        for(const auto &dist_pair : first_dists[0].first)
        {
            minimumOrder[min_index_count++] = dist_pair.second;
            // std::cout<< "( "<< dist_pair.first<< " : "<< dist_pair.second<< " )";
        }
        // std::cout<< std::endl;

        if(minimumOrder.size() != _sites.size())
        {
            throw("Error minimumorder.size != _sites.size()");
        }
        // for (const auto &dist_pair_pair : first_dists)
        // {
        //     minimumOrder[min_index_count++] = dist_pair_pair.second.second;
        // }



        setThisOrder(minimumOrder);
        if( _distances[0] != first_dists[0].first[0].first)
        {
            // std::cout<<"error distances not as expected, expected: "
            //          <<first_dists[0].first[0].first<< ". Gotten: "
            //          << _distances[0]<<std::endl;
        }
    }

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

        // std::cout<<"set this order before: "<<std::endl;
        // for(const auto tup : dist_indices)
        // {
        //     std::cout<<std::get<0>(tup)<< " "<<std::get<1>(tup)<< " "<<std::get<2>(tup)<<std::endl;
        // }

        std::sort(dist_indices.begin(), dist_indices.end());

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

    Sorted double vector
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
