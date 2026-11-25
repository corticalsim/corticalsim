#ifndef EDGE_HPP
#define EDGE_HPP

#include "types.hpp"
#include "linalg.hpp"

class Edge
{
  public:

    int exclude_point;

    double p_cat;
    double edg_ang_norm;
    double xy_rotation_angle;
    double edg_angle;
    double edg_bend_angle;

    Eigen::Vector2d mid_point;
    Eigen::Vector2d b;
    Eigen::Matrix2d a;
    vector<int> dir;
    Eigen::Vector2d dir2_d;
    vector<int> orientation;

    map<int, int> orientation_map;
    map<int, int> exclude_point_map;

    Edge()
    {
        exclude_point = 0;

        p_cat = 0.0;
        edg_ang_norm = 0.0;
        xy_rotation_angle = 0.0;
        edg_angle = 0.0;
        edg_bend_angle = 0.0;

        mid_point << 0.0, 0.0;

        a << 0.0, 0.0, 0.0, 0.0;

        b << 0.0, 0.0;

        orientation.push_back(0);
        orientation.push_back(0);
    }
};

#endif // EDGE_HPP
