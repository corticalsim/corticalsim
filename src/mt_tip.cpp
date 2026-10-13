#include "mt_tip.hpp"
#include "trajectory.hpp"
#include "microtubule.hpp"
#include "region.hpp"

MTTip::MTTip(Microtubule* m, TrajectoryVector tv, double v, DeterministicQueue* q, double queueV):
    mt(m),
    trajectory(tv.trajectory),
    dir(tv.dir),
    velocity(v),
    event(m, q, queueV),
    nextCollision(tv.trajectory->intersections.end())
{
// WARNING: Do not use the mt pointer in this function, as it has not been initialized
#ifdef DBG_MTS
    cout << "DBG/MTS: MTTip created.\n";
#endif

    return;
}

MTTip::~MTTip()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: MTTip destroyed.\n";
#endif

    return;
}

void MTTip::initialize()
{
    // register the tip for notification in the current trajectory
    notificationTag = trajectory->registerForNotifications(this);

    // get the region tag of the tip
    regionTag = trajectory->base.region->registerOnRegion(this, type(), mt->type);

    // locate the next intersection position of the tip
    locateIntersection();

    // determine the kind of event the tip is going to encounter
    determineEvent();

    return;
}

void MTTip::unlinkFromTrajectory()
{
    // unregister the tip from current region
    trajectory->base.region->unregisterFromRegion(regionTag, type(), mt->type);

    // also unregister its notification from current trajectory
    trajectory->unregisterForNotifications(notificationTag);
    return;
}

void MTTip::locateIntersection()
{
// go through the intersection list looking for the first event that can take place
#ifdef DBG_MTS
    cout << "DBG/MTS: MTTip ::locateIntersection()\n";
#endif

    // find the first intersection after current 'position()'
    nextCollision = trajectory->intersections.upper_bound(position());

    IntersectionItr temp = nextCollision;

    // if necessary, step down to the intersection on the other side of the tip [minus end can have velocity zero (=0)]
    if (((static_cast<double>(dir) * velocity) < 0) || ((velocity == 0) && (dir == ::forward)))
    {
        nextCollision--;
    }

    // if direction of the tip is backward then only step down the segment interator
    if (dir == backward)
    {
        temp--;
    }

    // plus tip (must be on the end of top segment)
    if (type() == t_plus)
    {
        mt->segments.last()->endItr = temp;
    }

    // minus tip (must be on the start of bottom segment
    else
    {
        mt->segments.first()->startItr = temp;
    }

    return;
}

void MTTip::advanceIntersection()
{
    // next collission site in front
    if ((static_cast<double>(dir) * velocity > 0) || ((velocity == 0) && (dir == backward)))
    {
#ifdef DBG_ASSERT
        if (nextCollision == trajectory->wallEnd())
        {
            cerr << "DBG/ASSERT: ERROR: Cannot advance position of tip that is at the boundary.\n";
            return;
        }
#endif

        nextCollision++;
    }

    // next collission site in back
    else
    {
#ifdef DBG_ASSERT
        if (nextCollision == trajectory->wallBegin())
        {
            cerr << "DBG/ASSERT: ERROR: Cannot advance position of tip that is at the boundary.\n";
            return;
        }
#endif

        nextCollision--;
    }

    IntersectionItr temp = nextCollision;

    // for minus end velocity (if direction of the tip is backward then step down the segment interator else step up)
    if (velocity <= 0)
    {
        if (dir == backward)
        {
            temp--;
        }
        else
        {
            temp++;
        }
    }

    // plus tip (must be on the end of top segment)
    if (type() == t_plus)
    {
        mt->segments.last()->endItr = temp;
    }

    // minus tip (must be on the start of bottom segment
    else
    {
        mt->segments.first()->startItr = temp;
    }

    return;
}

void MTTip::determineEvent()
{
    double eventPos(0.0);
    DeterministicEventType eventType;

    // get the current position of the tip on the trajectory
    double pos = position();

    // a tip that is polymerizing
    if (velocity > 0)
    {
        // next collision on wall begenning (at the bottom of the trajectory)
        if (nextCollision == trajectory->wallBegin())
        {
            eventType = ev_wall;
            eventPos = 0;
        }

        // next collision on wall end (at the top of the trajectory)
        else if (nextCollision == trajectory->wallEnd())
        {
            eventType = ev_wall;
            eventPos = trajectory->length;
        }

        // next collision, some where at the middle of trajectory
        else
        {
            eventType = ev_collision;
            eventPos = nextCollision->first;
        }
    }

    // a tip that is either paused or depolymerizing
    else
    {
        // next collision on wall begenning (at the bottom of the trajectory)
        if (nextCollision == trajectory->wallBegin())
        {
            // the segment will be removed
            eventType = ev_end_of_segment;
            eventPos = 0;
        }

        // next collision on wall end (at the top of the trajectory)
        else if (nextCollision == trajectory->wallEnd())
        {
            // the segment will be removed
            eventType = ev_end_of_segment;
            eventPos = trajectory->length;
        }

        // plus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->plus))
                 && (nextCollision == mt->segments.last()->startItr))
        {
            // the segment will be removed
            eventType = ev_end_of_segment;
            eventPos = mt->segments.last()->startItr->first;
        }

        // minus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->minus))
                 && (nextCollision == mt->segments.first()->endItr))
        {
            // the segment will be removed
            eventType = ev_end_of_segment;
            eventPos = mt->segments.first()->endItr->first;
        }

        // next collision is some where at the middle of trajectory and a back-track event will occur
        else
        {
            eventType = ev_backtrack;
            eventPos = nextCollision->first;
        }
    }

    // store the event with: (a) type of event and (b) distance to travel
    event.pushOnQueue((eventPos - pos) * static_cast<double>(dir), eventType);

    // transfer the value of the calculated event position to the next event position
    nextEventPos = eventPos;

#ifdef DBG_EVENT
    if (event.type != ev_none)
    {
        cout << "DBG/EVENT: event created. type: " << eventType << ", distance: " << (eventPos - pos) * dir << "\n";
        cout << "tip position: " << position() << ", tip velocity: " << dir * velocity
             << ", trajectory length: " << trajectory->length << "\n";
        if (nextCollision != trajectory->intersections.end())
        {
            cout << "event position: " << nextCollision->first << "\n";
        }
    }
#endif

    return;
}

void MTTip::notifyRemove(IntersectionItr& oldIs)
{
    // update the length of MT
    mt->updateLength();

    // scheduled next collision is at this invalid intersection point, so skip it to reach the next intersection point
    if (oldIs == nextCollision)
    {
        advanceIntersection();
        determineEvent();
    }

    // get the tip start/end iterator
    IntersectionItr& ref = (type() == t_plus) ? mt->segments.last()->endItr : mt->segments.first()->startItr;

    // link the tip start/end iterator to the new inserted iterator
    if (ref == oldIs)
    {
        if (dir == ::forward)
        {
            ref++;
        }
        else
        {
            ref--;
        }
    }

    return;
}

void MTTip::notifyInsert(IntersectionItr& newIs)
{

    // update the MT
    mt->updateLength();

    IntersectionItr temp = newIs;

    // next collission site in front
    if ((static_cast<double>(dir) * velocity > 0) || ((velocity == 0) && (dir == backward)))
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new
        // intersection
        if ((++temp == nextCollision) && (newIs->first > position()))
        {
            nextCollision = newIs;
            determineEvent();
        }
    }

    // next collission site in back
    else
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new
        // intersection
        if ((--temp == nextCollision) && (newIs->first < position()))
        {
            nextCollision = newIs;
            determineEvent();
        }
    }

    temp = newIs;

    // get the tip start/end iterator
    IntersectionItr& refItr = (type() == t_plus) ? mt->segments.last()->endItr : mt->segments.first()->startItr;

    // get the tip start/end position
    double& refPos = (type() == t_plus) ? mt->segments.last()->end : mt->segments.first()->start;

    // link the tip start/end iterator to the new inserted iterator
    if (dir == ::forward)
    {
        if ((++temp == refItr) && (newIs->first > refPos))
        {
            refItr = newIs;
        }
    }

    else
    {
        if ((--temp == refItr) && (newIs->first < refPos))
        {
            refItr = newIs;
        }
    }

    return;
}

TipType MTTip::type()
{
    // get tip type of MT
    if (this == &(mt->plus))
    {
        return t_plus;
    }
    else
    {
        return t_minus;
    }
}

double MTTip::position()
{
    // get current tip (top) position on trajectory
    if (this == &(mt->minus))
    {
        return mt->segments.first()->start;
    }
    else
    {
        return mt->segments.last()->end;
    }
}

double MTTip::otherPosition()
{
    // get the (bottom) position of an active segment
    if (this == &(mt->minus))
    {
        return mt->segments.first()->end;
    }
    else
    {
        return mt->segments.last()->start;
    }
}

Segment& MTTip::segment()
{
    // get the current segement of the tip
    if (this == &(mt->minus))
    {
        return *(mt->segments.first());
    }
    else
    {
        return *(mt->segments.last());
    }
}

void MTTip::switchTrajectory(Trajectory* newTr, Direction d, IntersectionItr intersect, bool advance)
{

    // copy old trajectory of the tip (required to unregister the tip from old trajectory region)
    Trajectory* oldTr(trajectory);

    // assign new trajectory to the tip
    trajectory = newTr;

    // assign new direction to the tip
    dir = d;
    nextCollision = intersect;

    // set the tip at the next intersection point
    if (advance)
    {
        advanceIntersection();
    }

    // assign an event to the tip
    determineEvent();

    // unregister the tip from old trajectory
    oldTr->base.region->unregisterFromRegion(regionTag, type(), mt->type);

    // find the new regiontag of the tip
    regionTag = newTr->base.region->registerOnRegion(this, type(), mt->type);

    // copy old notification tag of the tip (required to unregister the tip from old trajectory notification list)
    TrjMTTipTag tempTag = notificationTag;

    // assign new notification tag to the tip
    notificationTag = newTr->registerForNotifications(this);

    // unregister notification of the tip from old trajectory
    oldTr->unregisterForNotifications(tempTag);

    return;
}
