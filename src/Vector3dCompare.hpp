#pragma once
#include <eigen3/Eigen/Dense>
using namespace Eigen;

struct Vector3dCompare
{
    bool operator()(const Vector3d &lhs, const Vector3d &rhs) const
    {
        for (int i = 0; i < 3; i++)
        {
            if (lhs[i] == rhs[i])
            {
                continue;
            }
            else
            {
                return lhs[i] < rhs[i];
            }
        }
        //all equal return false
        return false;
    }
};

struct NeighborPairCompare
{
    bool operator()(const std::pair<int, Vector3d> &lhs, const std::pair<int, Vector3d> &rhs) const
    {
        if (lhs.first < rhs.first)
        {
            return true;
        }
        if (lhs.first > rhs.first)
        {
            return false;
        }

        for (int i = 0; i < 3; i++)
        {
            if (lhs.second[i] == rhs.second[i])
            {
                continue;
            }
            else
            {
                return lhs.second[i] < rhs.second[i];
            }
        }
        //all equal return false
        return false;
    }
};