#pragma once
#include <vector>
#include <boost/functional/hash.hpp>
using boost::hash;
using boost::hash_combine;
using boost::hash_value;



struct VectorHash {
    size_t operator()(const std::vector<int>& v) const 
    {
        size_t seed = 0;
        for (const int &i : v) {
            hash_combine(seed, hash_value(i));
        }
        return seed;
    }
};