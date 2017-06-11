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
        std::vector<std::vector<double>> first_dists;

        for(int i  = 0; i < _sites.size(); i++)
        {
            first_dists.push_back(getDistsToSite(i));
        }
        //ordering


    }

    /*
    Get all dists that origin from site i

    Sorted double vector
    */
    std::vector<double> getDistsToSite(int i)
    {
        std::vector<double> dists;
        int counter = 0;
        for(int k; k < _sites.size(); k++)
        {
            for(int l=k; l < _sites.size(); l++)
            {
                if (k== i or l == i)
                {
                    dists.push_back(_distances[counter]);
                }
                counter++;                
            }
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

        std::vector<std::tuple<int,int, double>> dist_indices;
        int counter = 0;
        for(int k; k < _sites.size(); k++)
        {
            for(int l=k; l < _sites.size(); l++)
            {
                int first = k;
                int second = l;
                if( first == i){ first = j; }
                if( first == j){ first = i; }
                if( second == i){ second = j; }
                if( second == j){ second = i; }
                
                if( second < first)
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
        for(const auto &tup : dist_indices)
        {
            _distances[counter++] = std::get<2>(tup);
        }
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

  private:
    std::vector<int> _sites;
    std::vector<double> _distances;
    std::map<std::vector<int>, int> _element_counts;
};