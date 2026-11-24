#ifndef INTERSECTION_HPP
#define INTERSECTION_HPP

#include "types.hpp"

class Intersection
// this object describes the intersection of one trajectory with another, from
// one side. since it is created twice for every intersection, it is by far the
// most memory-intensive object and should be kept as small as possible
{
  public:

    int m_occupancy;
    Trajectory* m_other_trajectory;
    IntersectionItr m_mirror;
#ifdef CROSS_SEV
    OccupiedIntersection* m_occupied_list_ptr;
#endif
};

class PointATedge
{
  public:

    double m_x;
    double m_y;
    double m_z;

    int m_next_element;

    PointATedge(double x1, double y1, double z1, int n1)
    {
        m_x = x1;
        m_y = y1;
        m_z = z1;
        m_next_element = n1;
    }
};

#ifdef CROSS_SEV
class OccupiedIntersection: public CompactListItem<OccupiedIntersection>
{
  public:

    OccupiedIntersection(IntersectionItr is):
        m_intersection_to_cut(is) {};
    IntersectionItr m_intersection_to_cut;
};
#endif

#endif // INTERSECTION_HPP
