#include "trajectory.hpp"
#include "region.hpp"
#include "microtubule.hpp"
#include "segment.hpp"
#include "system.hpp"

Trajectory::Trajectory(SurfaceVector base_vec, vector<PointATedge> end_p, double l):
    m_base(base_vec),
    m_length(l),
    m_prev_tr(0, ::FORWARD, nullptr),
    m_next_tr(0, ::FORWARD, nullptr),
    m_prev_tr_cos_angle(1.0),
    m_next_tr_cos_angle(1.0)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory created\n";
#endif

    // assign end points of the trajectory
    for (int i = 0; i < 2; i++)
    {
        m_end_point.push_back(end_p[i]);
    }

    // increase total number of trajectories in the geometry
    m_base.m_region->m_geometry->m_system->m_count_trajectories++;

    // initialize the intersection list of the trajectory
    m_intersections.insert(pair<double, Intersection>(-1.0, Intersection()));

    return;
}

Trajectory::~Trajectory()
{
    // decrease total number of trajectories in the geometry
    m_base.m_region->m_geometry->m_system->m_count_trajectories--;
    return;
}

bool Trajectory::integrity_check()
{
    bool valid(true);

    if (m_prev_tr.m_trajectory != nullptr)
    {
        if (m_prev_tr.m_dir == BACKWARD)
        {
            if (m_prev_tr.m_trajectory->m_next_tr.m_trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
        else
        {
            if (m_prev_tr.m_trajectory->m_prev_tr.m_trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
    }
    if (m_next_tr.m_trajectory != nullptr)
    {
        if (m_next_tr.m_dir == ::FORWARD)
        {
            if (m_next_tr.m_trajectory->m_prev_tr.m_trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
        else
        {
            if (m_next_tr.m_trajectory->m_next_tr.m_trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
    }

    IntersectionItr is(wall_begin());

    for (size_t is_count = 1; is_count < m_intersections.size(); is_count++)
    {
        is++;
        if (is->second.m_mirror->second.m_mirror != is)
        {
            cout << "Integrity: intersection mirroring broken.\n";
            valid = false;
        }
    }

    list<Segment*>::iterator seg_itr(m_segments.begin());

    while (seg_itr != m_segments.end())
    {
        is = (**seg_itr).m_start_itr;
        if ((**seg_itr).m_dir == ::FORWARD)
        {
            while (++is != (**seg_itr).m_end_itr)
            {
                is->second.m_occupancy--;
            }
        }
        else
        {
            while (--is != (**seg_itr).m_end_itr)
            {
                is->second.m_occupancy--;
            }
        }
        seg_itr++;
    }
    // now, check whether all occupancy counts are zero
    is = m_intersections.begin();
    while (++is != m_intersections.end())
    {
        if (is->second.m_occupancy != 0)
        {
            cerr << "Detected occupancy " << is->second.m_occupancy << " (should be zero).\n";
            valid = false;
        }
    }

    seg_itr = m_segments.begin();
    while (seg_itr != m_segments.end())
    {
        is = (**seg_itr).m_start_itr;
        if ((**seg_itr).m_dir == ::FORWARD)
        {
            while (++is != (**seg_itr).m_end_itr)
            {
                is->second.m_occupancy++;
            }
        }
        else
        {
            while (--is != (**seg_itr).m_end_itr)
            {
                is->second.m_occupancy++;
            }
        }
        seg_itr++;
    }

    if (!valid)
    {
        cerr << "Invalid occupancy count found during integrity check.\n";
    }

    return valid;
}

double Trajectory::segment_length()
{
    double sum_length(0.0);
    list<Segment*>::iterator seg(m_segments.begin());

    // calculate total length of all the segments present in the trajectory
    while (seg != m_segments.end())
    {
        sum_length += (**seg).length();
        seg++;
    }
    return sum_length;
}

double Trajectory::covered_length()
{
    list<Segment*>::iterator seg(m_segments.begin());
    list<double> start_points;
    list<double> end_points;

    // store the start and end point locations of segments in two seperate list
    while (seg != m_segments.end())
    {
        start_points.push_back(min((**seg).m_start, (**seg).m_end));
        end_points.push_back(max((**seg).m_start, (**seg).m_end));
        seg++;
    }

    // make sure the last value of start-point-list is maximum of all values from both lists
    start_points.push_back(VERY_LARGE);

    // sort the distances
    start_points.sort();
    end_points.sort();

    double previous_pos(0.0);
    double current_pos(0.0);
    int occupancy(0);
    double c_length(0.0);
    list<double>::iterator stepup(start_points.begin());
    list<double>::iterator stepdown(end_points.begin());

    while (stepdown != end_points.end())
    {
        current_pos = min(*stepup, *stepdown);

        // sum up only once and if the trajectory is covered, i.e. avoid multiple addition for bundles
        if (occupancy >= 1)
        {
            c_length += current_pos - previous_pos;
        }

        previous_pos = current_pos;

        // occupancy increasing, chances of the presence of bundle increasing
        if ((*stepup) < (*stepdown))
        {
            stepup++;
            occupancy++;
        }

        // occupancy decreasing, chances of the presence of bundle decresing
        else
        {
            stepdown++;
            occupancy--;
        }
    }

    return c_length;
}

void Trajectory::invalidate_intersection(IntersectionItr& old_is)
{
    list<MTTip*>::iterator tip(m_notification_list.begin());

    // notify all tips on the current trajectory that a previous insertion is invalidated
    while (tip != m_notification_list.end())
    {
        (**tip).notify_remove(old_is);
        tip++;
    }

    return;
}

void Trajectory::new_intersection(IntersectionItr& new_is)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory::newIntersection() called.\n";
#endif

    list<Segment*>::iterator seg(m_segments.begin());

    int occupancy(0);

    // find whether the new intersection is a cross-intersection for the segment
    while (seg != m_segments.end())
    {
        (**seg).m_mt->update_length();

        if ((**seg).crosses_intersection(new_is))
        {
            occupancy++;
        }

        seg++;
    }

    // assign complemetary occupancy for the new intersection
    new_is->second.m_occupancy = occupancy;

    // notify all tips on the current trajectory that a new insertion has happened
    list<MTTip*>::iterator tip(m_notification_list.begin());
    while (tip != m_notification_list.end())
    {
        (**tip).notify_insert(new_is);
        tip++;
    }

    return;
}

int Trajectory::difference_sign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2)
{
    // returns the sign of 'pos1 - pos2'. If the two positions are equal, the iterators are used instead to determine
    // order.
    if (pos1 - pos2 > ZERO_CUTOFF)
    {
        return 1;
    }

    else if (pos1 - pos2 < -ZERO_CUTOFF)
    {
        return -1;
    }

    else
    {
        if (itr1 == itr2)
        {
            return 0;
        }
        if (itr1 == wall_end())
        {
            return 1;
        }

        do
        {
            itr2++;
        } while ((itr2 != wall_end()) && (itr1 != itr2));

        if (itr1 == itr2)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
}

TrjSegmentTag Trajectory::insert_segment(Segment* s)
{
#ifdef DBG_ASSERT
    // check whether the size of the inserted segment equals zero, for insertion of non-zero segments, take care to
    // update the occupancy numbers
    if (abs(s->end - s->start) > ZERO_CUTOFF)
    {
        cerr << "DBG/ASSERT: ERROR: Not permitted to insert segment with non-zero length.\n";
        exit(-1);
    }
#endif

    // insert a new segment to the trajectory
    return m_segments.insert(m_segments.end(), s);
}

void Trajectory::remove_segment(TrjSegmentTag s)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory::removeSegment() called.\n";
#endif

    // remove segment
    m_segments.erase(s);

    // if the associated trajectory has no segement and no tips, then remove the whole trajectory
    conditional_remove();

    return;
}

TrjMTTipTag Trajectory::register_for_notifications(MTTip* p_tip)
{
    // a new tip arrived on the trajectory, add it to the notification list
    return m_notification_list.insert(m_notification_list.end(), p_tip);
}

void Trajectory::unregister_for_notifications(TrjMTTipTag tag)
{
    // a tip is leaving the trajectory, remove it from the notification list
    m_notification_list.erase(tag);

    // if the associated trajectory has no tips and no segments, then remove the whole trajectory
    conditional_remove();
    return;
}

TrajectoryVector Trajectory::next_trajectory(Direction dir)
{
    if (dir == ::FORWARD)
    {
        // direction forward and a next trajectory is waiting to take the tip on board
        if (m_next_tr.m_trajectory != nullptr)
        {
#ifdef DBG_GEOMETRY
            cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from "
                 << RegionTypeText[base.region->type] << " to " << RegionTypeText[nextTr.trajectory->base.region->type]
                 << "]\n";
#endif

            return m_next_tr;
        }
    }

    else
    {
        // direction backward and a previous trajectory is waiting to take the on board
        if (m_prev_tr.m_trajectory != nullptr)
        {
#ifdef DBG_GEOMETRY
            cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from "
                 << RegionTypeText[base.region->type] << " to " << RegionTypeText[prevTr.trajectory->base.region->type]
                 << "]\n";
#endif

            return m_prev_tr;
        }
    }

    // there is no next/previous trajectory to take the tip, so extend (create a new) trajectory
    return m_base.m_region->m_geometry->extend_trajectory(this, dir);
}

void Trajectory::conditional_remove()
{
    // if there is no segments and no tips, associted with a trajectory, then remove it rigth way
    if (m_segments.empty() && (m_notification_list.size() == 0))
    {
        m_base.m_region->remove_trajectory(this);
    }

    return;
}
