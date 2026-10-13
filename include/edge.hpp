#ifndef EDGE_HPP
#define EDGE_HPP

#include "types.hpp"

class Edge
{
  public:

    int excludePoint;

    double pCat;
    double edgAngNorm;
    double xyRotationAngle;
    double edgAngle;
    double edgBendAngle;

    Vector2d midPoint;
    Vector2d b;
    Matrix2d A;
    vector<int> dir;
    Vector2d dir2D;
    vector<int> orientation;

    map<int, int> orientationMap;
    map<int, int> excludePointMap;

    Edge()
    {
        excludePoint = 0;

        pCat = 0.0;
        edgAngNorm = 0.0;
        xyRotationAngle = 0.0;
        edgAngle = 0.0;
        edgBendAngle = 0.0;

        midPoint << 0.0, 0.0;

        A << 0.0, 0.0, 0.0, 0.0;

        b << 0.0, 0.0;

        orientation.push_back(0);
        orientation.push_back(0);
    }
};

#endif // EDGE_HPP
