#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include "types.hpp"
#include "trajectory.hpp"

class Segment: public DLBaseItem<Segment>
{
  public:

    Microtubule* mt;
    Trajectory* const trajectory;
    TrjSegmentTag trajectory_tag;
    double nucleation_time;

    // -+ direction with respect to trajectory orientation
    Direction dir;
    double start;
    double end;
    IntersectionItr start_itr;
    IntersectionItr end_itr;

    // constructs a segment as part of a microtubule at a specified vector
    // location
    Segment(Microtubule*, TrajectoryVector&);

    ~Segment();

    double length() { return abs(end - start); }

    bool is_last_in_mt();
    bool is_first_in_mt();
    bool crosses_intersection(IntersectionItr& is);
};

#endif // SEGMENT_HPP
