#ifndef EDGE_HPP
#define EDGE_HPP

#include "types.hpp"
#include "linalg.hpp"

class Edge
{
  public:

    int m_exclude_point;

    double m_p_cat;
    double m_edg_ang_norm;
    double m_xy_rotation_angle;
    double m_edg_angle;
    double m_edg_bend_angle;

    Eigen::Vector2d m_mid_point;
    Eigen::Vector2d m_b;
    Eigen::Matrix2d m_a;
    vector<int> m_dir;
    Eigen::Vector2d m_dir2_d;
    vector<int> m_orientation;

    map<int, int> m_orientation_map;
    map<int, int> m_exclude_point_map;

    Edge()
    {
        m_exclude_point = 0;

        m_p_cat = 0.0;
        m_edg_ang_norm = 0.0;
        m_xy_rotation_angle = 0.0;
        m_edg_angle = 0.0;
        m_edg_bend_angle = 0.0;

        m_mid_point << 0.0, 0.0;

        m_a << 0.0, 0.0, 0.0, 0.0;

        m_b << 0.0, 0.0;

        m_orientation.push_back(0);
        m_orientation.push_back(0);
    }
};

#endif // EDGE_HPP
