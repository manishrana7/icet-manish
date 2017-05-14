#pragma once
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <vector>

/**
Design approach:
    input pair neighbors and calculate higher order neighbors
    using set intersection.
*/

class ManybodyNeighborlist
{
    public:
    ManybodyNeighborlist(const <td::vector<double> cutoffs)
    {
        cutoffs = cutoffs;
    }
    void build();

    


    private:
    _cutoffs;

};