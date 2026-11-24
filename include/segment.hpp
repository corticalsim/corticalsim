#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include "types.hpp"
#include "trajectory.hpp"

class Segment: public DLBaseItem<Segment>
{
  public:

    Microtubule* m_mt;
    Trajectory* const m_trajectory;
    TrjSegmentTag m_trajectory_tag;
    double m_nucleation_time;

    // -+ direction with respect to trajectory orientation
    Direction m_dir;
    double m_start;
    double m_end;
    IntersectionItr m_start_itr;
    IntersectionItr m_end_itr;

    // constructs a segment as part of a microtubule at a specified vector
    // location
    Segment(Microtubule*, TrajectoryVector&);

    ~Segment();

    double length() { return abs(m_end - m_start); }

    bool is_last_in_mt();
    bool is_first_in_mt();
    bool crosses_intersection(IntersectionItr& is);
};

#endif // SEGMENT_HPP
