#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "types.hpp"
#include "linalg.hpp"

struct Vertices
{

    double m_x;
    double m_y;
    double m_z;

    Vertices(double x1, double y1, double z1)
    {
        m_x = x1;
        m_y = y1;
        m_z = z1;
    }
};

struct BendingOperator
{

    double m_cos_theta;
    double m_sin_theta;

    Eigen::Vector3d m_axis;

    BendingOperator()
    {
        m_cos_theta = 0.0;
        m_sin_theta = 0.0;

        m_axis << 0.0, 0.0, 0.0;
    }
};

struct SurfaceVector
{
    double m_x;
    double m_y;
    double m_z;
    double m_angle;
    double m_tv_pos;

    Region* m_region;

    SurfaceVector()
    {
        m_x = 0.0;
        m_y = 0.0;
        m_z = 0.0;

        m_angle = 0.0;
        m_tv_pos = 0.0;
    }
};

struct ElementList
{

    int m_mid_point;
    double m_length;
    vector<int> m_nodes;
    vector<int> m_shared_element;
    vector<int> m_shared_edge;
    double m_edge3_d_angle;
    Eigen::Vector3d m_intersection_by_plane;

    ElementList()
    {
        m_mid_point = 0;
        m_length = 0.0;
        m_edge3_d_angle = 0.0;
        m_intersection_by_plane << 0.0, 0.0, 0.0;

        for (int i = 0; i < 2; i++)
        {
            m_nodes.push_back(0);
            m_shared_element.push_back(0);
            m_shared_edge.push_back(0);
        }
    }
};

#endif // SURFACE_HPP
