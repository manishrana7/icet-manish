#pragma once

#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include "LatticeNeighbor.hpp"
#include "Structure.hpp"


namespace icet
{

/// Returns the geomtetrical radius from the vectors of latticeneighbors and a unitcell
double getGeometricalRadius(const std::vector<LatticeNeighbor> &, const Structure &);



}
