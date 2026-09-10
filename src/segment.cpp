#include "segment.hpp"
#include "system.hpp"
#include "microtubule.hpp"

Segment::Segment(Microtubule* m, TrajectoryVector& tv):
    mt(m),
    trajectory(tv.trajectory),
    trajectoryTag(tv.trajectory->segments.end()),
    nucleationTime(m->system->systemTime + m->system->systemTimeOffset),
    dir(tv.dir),
    start(tv.pos),
    end(tv.pos),
    startItr(tv.trajectory->wallEnd()),
    endItr(tv.trajectory->wallEnd())
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Segment created.\n";
    cout << "begin position: " << start << ", end position: " << end << "\n";
#endif

    // insert this segment to the associated MT trajectory
    trajectoryTag = tv.trajectory->insertSegment(this);

    // increase total number of segments for the associated MT
    mt->system->countSegments++;
    return;
}

Segment::~Segment()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Segment destroyed.\n";
#endif

    // decrease total number of segments for the associated MT
    mt->system->countSegments--;

    // if this is the last segment, remove the tip references...
    trajectory->removeSegment(trajectoryTag);
    return;
}

bool Segment::isLastInMT()
{
    // check whether this is the last segment of the MT
    return (this == mt->segments.last());
}

bool Segment::isFirstInMT()
{
    // check whether this is the first segment of the MT
    return (this == mt->segments.first());
}

bool Segment::crossesIntersection(IntersectionItr& is)
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
        int sign1 = trajectory->differenceSign(
        (this == mt->segments.last()) ? mt->plus.nextCollision : endItr, end, is, is->first);
        int sign2 = trajectory->differenceSign(
        (this == mt->segments.first()) ? mt->minus.nextCollision : startItr, start, is, is->first);

        if (sign1 * sign2 == -1)
        {
            return true;
        }
    }

    return false;
}
