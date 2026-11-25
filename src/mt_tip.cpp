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
    next_collision(tv.trajectory->intersections.end())
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
    notification_tag = trajectory->register_for_notifications(this);

    // get the region tag of the tip
    region_tag = trajectory->base.region->register_on_region(this, type(), mt->type);

    // locate the next intersection position of the tip
    locate_intersection();

    // determine the kind of event the tip is going to encounter
    determine_event();

    return;
}

void MTTip::unlink_from_trajectory()
{
    // unregister the tip from current region
    trajectory->base.region->unregister_from_region(region_tag, type(), mt->type);

    // also unregister its notification from current trajectory
    trajectory->unregister_for_notifications(notification_tag);
    return;
}

void MTTip::locate_intersection()
{
// go through the intersection list looking for the first event that can take place
#ifdef DBG_MTS
    cout << "DBG/MTS: MTTip ::locateIntersection()\n";
#endif

    // find the first intersection after current 'position()'
    next_collision = trajectory->intersections.upper_bound(position());

    IntersectionItr temp = next_collision;

    // if necessary, step down to the intersection on the other side of the tip [minus end can have velocity zero (=0)]
    if (((static_cast<double>(dir) * velocity) < 0) || ((velocity == 0) && (dir == Direction::FORWARD)))
    {
        next_collision--;
    }

    // if direction of the tip is backward then only step down the segment interator
    if (dir == Direction::BACKWARD)
    {
        temp--;
    }

    // plus tip (must be on the end of top segment)
    if (type() == TipType::PLUS)
    {
        mt->segments.last()->end_itr = temp;
    }

    // minus tip (must be on the start of bottom segment
    else
    {
        mt->segments.first()->start_itr = temp;
    }

    return;
}

void MTTip::advance_intersection()
{
    // next collission site in front
    if ((static_cast<double>(dir) * velocity > 0) || ((velocity == 0) && (dir == Direction::BACKWARD)))
    {
#ifdef DBG_ASSERT
        if (nextCollision == trajectory->wallEnd())
        {
            cerr << "DBG/ASSERT: ERROR: Cannot advance position of tip that is at the boundary.\n";
            return;
        }
#endif

        next_collision++;
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

        next_collision--;
    }

    IntersectionItr temp = next_collision;

    // for minus end velocity (if direction of the tip is backward then step down the segment interator else step up)
    if (velocity <= 0)
    {
        if (dir == Direction::BACKWARD)
        {
            temp--;
        }
        else
        {
            temp++;
        }
    }

    // plus tip (must be on the end of top segment)
    if (type() == TipType::PLUS)
    {
        mt->segments.last()->end_itr = temp;
    }

    // minus tip (must be on the start of bottom segment
    else
    {
        mt->segments.first()->start_itr = temp;
    }

    return;
}

void MTTip::determine_event()
{
    double event_pos(0.0);
    DeterministicEventType event_type;

    // get the current position of the tip on the trajectory
    double pos = position();

    // a tip that is polymerizing
    if (velocity > 0)
    {
        // next collision on wall begenning (at the bottom of the trajectory)
        if (next_collision == trajectory->wall_begin())
        {
            event_type = DeterministicEventType::WALL;
            event_pos = 0;
        }

        // next collision on wall end (at the top of the trajectory)
        else if (next_collision == trajectory->wall_end())
        {
            event_type = DeterministicEventType::WALL;
            event_pos = trajectory->length;
        }

        // next collision, some where at the middle of trajectory
        else
        {
            event_type = DeterministicEventType::COLLISION;
            event_pos = next_collision->first;
        }
    }

    // a tip that is either paused or depolymerizing
    else
    {
        // next collision on wall begenning (at the bottom of the trajectory)
        if (next_collision == trajectory->wall_begin())
        {
            // the segment will be removed
            event_type = DeterministicEventType::END_OF_SEGMENT;
            event_pos = 0;
        }

        // next collision on wall end (at the top of the trajectory)
        else if (next_collision == trajectory->wall_end())
        {
            // the segment will be removed
            event_type = DeterministicEventType::END_OF_SEGMENT;
            event_pos = trajectory->length;
        }

        // plus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->plus))
                 && (next_collision == mt->segments.last()->start_itr))
        {
            // the segment will be removed
            event_type = DeterministicEventType::END_OF_SEGMENT;
            event_pos = mt->segments.last()->start_itr->first;
        }

        // minus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->minus))
                 && (next_collision == mt->segments.first()->end_itr))
        {
            // the segment will be removed
            event_type = DeterministicEventType::END_OF_SEGMENT;
            event_pos = mt->segments.first()->end_itr->first;
        }

        // next collision is some where at the middle of trajectory and a back-track event will occur
        else
        {
            event_type = DeterministicEventType::BACKTRACK;
            event_pos = next_collision->first;
        }
    }

    // store the event with: (a) type of event and (b) distance to travel
    event.push_on_queue((event_pos - pos) * static_cast<double>(dir), event_type);

    // transfer the value of the calculated event position to the next event position
    next_event_pos = event_pos;

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

void MTTip::notify_remove(IntersectionItr& old_is)
{
    // update the length of MT
    mt->update_length();

    // scheduled next collision is at this invalid intersection point, so skip it to reach the next intersection point
    if (old_is == next_collision)
    {
        advance_intersection();
        determine_event();
    }

    // get the tip start/end iterator
    IntersectionItr& ref = (type() == TipType::PLUS) ? mt->segments.last()->end_itr : mt->segments.first()->start_itr;

    // link the tip start/end iterator to the new inserted iterator
    if (ref == old_is)
    {
        if (dir == Direction::FORWARD)
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

void MTTip::notify_insert(IntersectionItr& new_is)
{

    // update the MT
    mt->update_length();

    IntersectionItr temp = new_is;

    // next collission site in front
    if ((static_cast<double>(dir) * velocity > 0) || ((velocity == 0) && (dir == Direction::BACKWARD)))
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new
        // intersection
        if ((++temp == next_collision) && (new_is->first > position()))
        {
            next_collision = new_is;
            determine_event();
        }
    }

    // next collission site in back
    else
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new
        // intersection
        if ((--temp == next_collision) && (new_is->first < position()))
        {
            next_collision = new_is;
            determine_event();
        }
    }

    temp = new_is;

    // get the tip start/end iterator
    IntersectionItr& ref_itr = (type() == TipType::PLUS) ? mt->segments.last()->end_itr : mt->segments.first()->start_itr;

    // get the tip start/end position
    double& ref_pos = (type() == TipType::PLUS) ? mt->segments.last()->end : mt->segments.first()->start;

    // link the tip start/end iterator to the new inserted iterator
    if (dir == Direction::FORWARD)
    {
        if ((++temp == ref_itr) && (new_is->first > ref_pos))
        {
            ref_itr = new_is;
        }
    }

    else
    {
        if ((--temp == ref_itr) && (new_is->first < ref_pos))
        {
            ref_itr = new_is;
        }
    }

    return;
}

TipType MTTip::type()
{
    // get tip type of MT
    if (this == &(mt->plus))
    {
        return TipType::PLUS;
    }
    else
    {
        return TipType::MINUS;
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

double MTTip::other_position()
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

void MTTip::switch_trajectory(Trajectory* new_tr, Direction d, IntersectionItr intersect, bool advance)
{

    // copy old trajectory of the tip (required to unregister the tip from old trajectory region)
    Trajectory* old_tr(trajectory);

    // assign new trajectory to the tip
    trajectory = new_tr;

    // assign new direction to the tip
    dir = d;
    next_collision = intersect;

    // set the tip at the next intersection point
    if (advance)
    {
        advance_intersection();
    }

    // assign an event to the tip
    determine_event();

    // unregister the tip from old trajectory
    old_tr->base.region->unregister_from_region(region_tag, type(), mt->type);

    // find the new regiontag of the tip
    region_tag = new_tr->base.region->register_on_region(this, type(), mt->type);

    // copy old notification tag of the tip (required to unregister the tip from old trajectory notification list)
    TrjMTTipTag temp_tag = notification_tag;

    // assign new notification tag to the tip
    notification_tag = new_tr->register_for_notifications(this);

    // unregister notification of the tip from old trajectory
    old_tr->unregister_for_notifications(temp_tag);

    return;
}
