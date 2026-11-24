#include "segment.hpp"
#include "system.hpp"
#include "microtubule.hpp"

Segment::Segment(Microtubule* m, TrajectoryVector& tv):
    m_mt(m),
    m_trajectory(tv.m_trajectory),
    m_trajectory_tag(tv.m_trajectory->m_segments.end()),
    m_nucleation_time(m->m_system->m_system_time + m->m_system->m_system_time_offset),
    m_dir(tv.m_dir),
    m_start(tv.m_pos),
    m_end(tv.m_pos),
    m_start_itr(tv.m_trajectory->wall_end()),
    m_end_itr(tv.m_trajectory->wall_end())
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Segment created.\n";
    cout << "begin position: " << start << ", end position: " << end << "\n";
#endif

    // insert this segment to the associated MT trajectory
    m_trajectory_tag = tv.m_trajectory->insert_segment(this);

    // increase total number of segments for the associated MT
    m_mt->m_system->m_count_segments++;
    return;
}

Segment::~Segment()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Segment destroyed.\n";
#endif

    // decrease total number of segments for the associated MT
    m_mt->m_system->m_count_segments--;

    // if this is the last segment, remove the tip references...
    m_trajectory->remove_segment(m_trajectory_tag);
    return;
}

bool Segment::is_last_in_mt()
{
    // check whether this is the last segment of the MT
    return (this == m_mt->m_segments.last());
}

bool Segment::is_first_in_mt()
{
    // check whether this is the first segment of the MT
    return (this == m_mt->m_segments.first());
}

bool Segment::crosses_intersection(IntersectionItr& is)
{
    // check whether the segment has a cross intersection with a trajectory
    double temp(0.0);
    temp = (is->first - m_end) * (is->first - m_start);

    if (temp < -ZERO_CUTOFF)
    {
        return true;
    }

    else if (temp < ZERO_CUTOFF)
    {
        int sign1 = m_trajectory->difference_sign(
        (this == m_mt->m_segments.last()) ? m_mt->m_plus.m_next_collision : m_end_itr, m_end, is, is->first);
        int sign2 = m_trajectory->difference_sign(
        (this == m_mt->m_segments.first()) ? m_mt->m_minus.m_next_collision : m_start_itr, m_start, is, is->first);

        if (sign1 * sign2 == -1)
        {
            return true;
        }
    }

    return false;
}
