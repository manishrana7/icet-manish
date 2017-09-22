#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <string>
#include "Structure.hpp"
#include "OrbitList.hpp"

using namespace Eigen;


/**
This is the cluster space object. 

It will have the definition of the cluster space a cluster expansion is based on.

*/

class ClusterSpace
{
public:
ClusterSpace(){};


///Generate the clustervector on the input structure
std::vector<double> generateClustervector(const Structure &) const;

private:
///Primitive cell/structure    
Structure _primitive_structure;

///Primitive orbitlist based on the structure and the global cutoffs
OrbitList _primitive_orbitlist;

///Unique id for this clusterspace
int clusterSpace_ID;

///The radial cutoffs for neigbhor inclusion. First element in vector is pairs, then triplets etc.
std::vector<double> _clusterCutoffs;

/**
This maps a orbit to other orbits i.e. of _equalSitesList[0] = [1,2,3] then orbit zero should be seen as equal to orbit 1,2 and three.
This will mean that when retrieving the cluster vector the zeroth element will be a combination of orbit 0,1,2 and 3.
*/
std::map<int,std::vector<int>> _equalSitesList;

//std::map<std::pair<int,int>, double> _defaultClusterFunction; 

///return the default clusterfunction of the input element when this site can have Mi elements
double _defaultClusterFunction(const int Mi, const int element);
};