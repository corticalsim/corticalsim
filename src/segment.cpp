#include "segment.hpp"
#include "system.hpp"
#include "microtubule.hpp"

Segment::Segment(Microtubule* m, TrajectoryVector& tv):
    mt(m),
    trajectory(tv.trajectory),
    trajectory_tag(tv.trajectory->segments.end()),
    nucleation_time(m->system->system_time + m->system->system_time_offset),
    dir(tv.dir),
    start(tv.pos),
    end(tv.pos),
    start_itr(tv.trajectory->wall_end()),
    end_itr(tv.trajectory->wall_end())
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Segment created.\n";
    cout << "begin position: " << start << ", end position: " << end << "\n";
#endif

    // insert this segment to the associated MT trajectory
    trajectory_tag = tv.trajectory->insert_segment(this);

    // increase total number of segments for the associated MT
    mt->system->count_segments++;
    return;
}

Segment::~Segment()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Segment destroyed.\n";
#endif

    // decrease total number of segments for the associated MT
    mt->system->count_segments--;

    // if this is the last segment, remove the tip references...
    trajectory->remove_segment(trajectory_tag);
    return;
}

bool Segment::is_last_in_mt()
{
    // check whether this is the last segment of the MT
    return (this == mt->segments.last());
}

bool Segment::is_first_in_mt()
{
    // check whether this is the first segment of the MT
    return (this == mt->segments.first());
}

bool Segment::crosses_intersection(IntersectionItr& is)
{
    // check whether the segment has a cross intersection with a trajectory
    double temp(0.0);
    temp = (is->first - end) * (is->first - start);

    if (temp < -ZERO_CUTOFF)
    {
        return true;
    }

    else if (temp < ZERO_CUTOFF)
    {
        int sign1 = trajectory->difference_sign(
        (this == mt->segments.last()) ? mt->plus.next_collision : end_itr, end, is, is->first);
        int sign2 = trajectory->difference_sign(
        (this == mt->segments.first()) ? mt->minus.next_collision : start_itr, start, is, is->first);

        if (sign1 * sign2 == -1)
        {
            return true;
        }
    }

    return false;
}
