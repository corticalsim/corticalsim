#include "mt_tip.hpp"
#include "trajectory.hpp"
#include "microtubule.hpp"
#include "region.hpp"

MTTip::MTTip(Microtubule* m, TrajectoryVector tv, double v, DeterministicQueue* q, double queue_v):
    mt(m),
    trajectory(tv.trajectory),
    dir(tv.dir),
    velocity(v),
    event(m, q, queue_v),
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
    if (((static_cast<double>(dir) * velocity) < 0) || ((velocity == 0) && (dir == ::FORWARD)))
    {
        nextCollision--;
    }

    // if direction of the tip is backward then only step down the segment interator
    if (dir == BACKWARD)
    {
        temp--;
    }

    // plus tip (must be on the end of top segment)
    if (type() == T_PLUS)
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
    if ((static_cast<double>(dir) * velocity > 0) || ((velocity == 0) && (dir == BACKWARD)))
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
        if (dir == BACKWARD)
        {
            temp--;
        }
        else
        {
            temp++;
        }
    }

    // plus tip (must be on the end of top segment)
    if (type() == T_PLUS)
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
            eventType = EV_WALL;
            eventPos = 0;
        }

        // next collision on wall end (at the top of the trajectory)
        else if (nextCollision == trajectory->wallEnd())
        {
            eventType = EV_WALL;
            eventPos = trajectory->length;
        }

        // next collision, some where at the middle of trajectory
        else
        {
            eventType = EV_COLLISION;
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
            eventType = EV_END_OF_SEGMENT;
            eventPos = 0;
        }

        // next collision on wall end (at the top of the trajectory)
        else if (nextCollision == trajectory->wallEnd())
        {
            // the segment will be removed
            eventType = EV_END_OF_SEGMENT;
            eventPos = trajectory->length;
        }

        // plus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->plus))
                 && (nextCollision == mt->segments.last()->startItr))
        {
            // the segment will be removed
            eventType = EV_END_OF_SEGMENT;
            eventPos = mt->segments.last()->startItr->first;
        }

        // minus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->minus))
                 && (nextCollision == mt->segments.first()->endItr))
        {
            // the segment will be removed
            eventType = EV_END_OF_SEGMENT;
            eventPos = mt->segments.first()->endItr->first;
        }

        // next collision is some where at the middle of trajectory and a back-track event will occur
        else
        {
            eventType = EV_BACKTRACK;
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

void MTTip::notifyRemove(IntersectionItr& old_is)
{
    // update the length of MT
    mt->updateLength();

    // scheduled next collision is at this invalid intersection point, so skip it to reach the next intersection point
    if (old_is == nextCollision)
    {
        advanceIntersection();
        determineEvent();
    }

    // get the tip start/end iterator
    IntersectionItr& ref = (type() == T_PLUS) ? mt->segments.last()->endItr : mt->segments.first()->startItr;

    // link the tip start/end iterator to the new inserted iterator
    if (ref == old_is)
    {
        if (dir == ::FORWARD)
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

void MTTip::notifyInsert(IntersectionItr& new_is)
{

    // update the MT
    mt->updateLength();

    IntersectionItr temp = new_is;

    // next collission site in front
    if ((static_cast<double>(dir) * velocity > 0) || ((velocity == 0) && (dir == BACKWARD)))
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new
        // intersection
        if ((++temp == nextCollision) && (new_is->first > position()))
        {
            nextCollision = new_is;
            determineEvent();
        }
    }

    // next collission site in back
    else
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new
        // intersection
        if ((--temp == nextCollision) && (new_is->first < position()))
        {
            nextCollision = new_is;
            determineEvent();
        }
    }

    temp = new_is;

    // get the tip start/end iterator
    IntersectionItr& refItr = (type() == T_PLUS) ? mt->segments.last()->endItr : mt->segments.first()->startItr;

    // get the tip start/end position
    double& refPos = (type() == T_PLUS) ? mt->segments.last()->end : mt->segments.first()->start;

    // link the tip start/end iterator to the new inserted iterator
    if (dir == ::FORWARD)
    {
        if ((++temp == refItr) && (new_is->first > refPos))
        {
            refItr = new_is;
        }
    }

    else
    {
        if ((--temp == refItr) && (new_is->first < refPos))
        {
            refItr = new_is;
        }
    }

    return;
}

TipType MTTip::type()
{
    // get tip type of MT
    if (this == &(mt->plus))
    {
        return T_PLUS;
    }
    else
    {
        return T_MINUS;
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

void MTTip::switchTrajectory(Trajectory* new_tr, Direction d, IntersectionItr intersect, bool advance)
{

    // copy old trajectory of the tip (required to unregister the tip from old trajectory region)
    Trajectory* oldTr(trajectory);

    // assign new trajectory to the tip
    trajectory = new_tr;

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
    regionTag = new_tr->base.region->registerOnRegion(this, type(), mt->type);

    // copy old notification tag of the tip (required to unregister the tip from old trajectory notification list)
    TrjMTTipTag tempTag = notificationTag;

    // assign new notification tag to the tip
    notificationTag = new_tr->registerForNotifications(this);

    // unregister notification of the tip from old trajectory
    oldTr->unregisterForNotifications(tempTag);

    return;
}
