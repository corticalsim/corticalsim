#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "types.hpp"

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

struct bendingOperator
{

    double cosTheta;
    double sinTheta;

    Vector3d axis;

    bendingOperator()
    {
        cosTheta = 0.0;
        sinTheta = 0.0;

        axis << 0.0, 0.0, 0.0;
    }
};

struct SurfaceVector
{
    double x;
    double y;
    double z;
    double angle;
    double tvPos;

    Region* region;

    SurfaceVector()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;

        angle = 0.0;
        tvPos = 0.0;
    }
};

struct elementList
{

    int midPoint;
    double length;
    vector<int> nodes;
    vector<int> sharedElement;
    vector<int> sharedEdge;
    double edge3DAngle;
    Vector3d intersectionByPlane;

    elementList()
    {
        midPoint = 0;
        length = 0.0;
        edge3DAngle = 0.0;
        intersectionByPlane << 0.0, 0.0, 0.0;

        for (int i = 0; i < 2; i++)
        {
            nodes.push_back(0);
            sharedElement.push_back(0);
            sharedEdge.push_back(0);
        }
    }
};

#endif // SURFACE_HPP
