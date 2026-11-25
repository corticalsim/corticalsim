#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "types.hpp"
#include "linalg.hpp"

struct Vertices
{

    double x;
    double y;
    double z;

    Vertices(double x1, double y1, double z1)
    {
        x = x1;
        y = y1;
        z = z1;
    }
};

struct BendingOperator
{

    double cos_theta;
    double sin_theta;

    Eigen::Vector3d axis;

    BendingOperator()
    {
        cos_theta = 0.0;
        sin_theta = 0.0;

        axis << 0.0, 0.0, 0.0;
    }
};

struct SurfaceVector
{
    double x;
    double y;
    double z;
    double angle;
    double tv_pos;

    Region* region;

    SurfaceVector()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;

        angle = 0.0;
        tv_pos = 0.0;
    }
};

struct ElementList
{

    int mid_point;
    double length;
    vector<int> nodes;
    vector<int> shared_element;
    vector<int> shared_edge;
    double edge3_d_angle;
    Eigen::Vector3d intersection_by_plane;

    ElementList()
    {
        mid_point = 0;
        length = 0.0;
        edge3_d_angle = 0.0;
        intersection_by_plane << 0.0, 0.0, 0.0;

        for (int i = 0; i < 2; i++)
        {
            nodes.push_back(0);
            shared_element.push_back(0);
            shared_edge.push_back(0);
        }
    }
};

#endif // SURFACE_HPP
