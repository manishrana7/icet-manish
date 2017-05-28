#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
using namespace Eigen;

namespace py = pybind11;

class ClusterCounts
{
    public:
    ClusterCounts();

    void count(const Structure & , XXX indices );


}