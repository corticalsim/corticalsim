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
    Trajectory* other_trajectory;
    IntersectionItr mirror;
#ifdef CROSS_SEV
    OccupiedIntersection* occupied_list_ptr;
#endif
};

class PointAtEdge
{
  public:

    double x;
    double y;
    double z;

    int next_element;

    PointAtEdge(double x1, double y1, double z1, int n1)
    {
        x = x1;
        y = y1;
        z = z1;
        next_element = n1;
    }
};

#ifdef CROSS_SEV
class OccupiedIntersection: public CompactListItem<OccupiedIntersection>
{
  public:

    OccupiedIntersection(IntersectionItr is):
        intersection_to_cut(is) {};
    IntersectionItr intersection_to_cut;
};
#endif

#endif // INTERSECTION_HPP
