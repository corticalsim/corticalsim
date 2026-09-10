#ifndef INTERSECTION_HPP
#define INTERSECTION_HPP

#include "types.hpp"

class Intersection
// this object describes the intersection of one trajectory with another, from
// one side. since it is created twice for every intersection, it is by far the
// most memory-intensive object and should be kept as small as possible
{
  public:

    int occupancy;
    Trajectory* otherTrajectory;
    IntersectionItr mirror;
#ifdef CROSS_SEV
    OccupiedIntersection* occupiedListPtr;
#endif
};

class PointATedge
{
  public:

    double x;
    double y;
    double z;

    int nextElement;

    PointATedge(double x1, double y1, double z1, int n1)
    {
        x = x1;
        y = y1;
        z = z1;
        nextElement = n1;
    }
};

#ifdef CROSS_SEV
class OccupiedIntersection: public CompactListItem<OccupiedIntersection>
{
  public:

    OccupiedIntersection(IntersectionItr is):
        intersectionToCut(is) {};
    IntersectionItr intersectionToCut;
};
#endif

#endif // INTERSECTION_HPP
