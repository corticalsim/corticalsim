#include "region.hpp"
#include "geometry.hpp"
#include "system.hpp"

TrajectoryVector Region::insert_trajectory(const SurfaceVector& s_vec)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Region::insertTrajectory() called\n";
#endif

    TrajectoryVector t_vec;
    double t_length(0.);
    SurfaceVector new_base(s_vec);

    vector<PointAtEdge> end_point;

    // gather coordinates of the trajectory to be created
    get_trajectory_coordinates(new_base, t_length, end_point, t_vec);

    // create the trajectory
    t_vec.trajectory = trajectories.create(new_base, end_point, t_length);

    // make a list of intersections of the created trajectory with all other trajectories
    make_intersection_list(t_vec.trajectory);

    // increase total number of intersections accordingly
    geometry->system->count_intersections += 2 * (t_vec.trajectory->intersections.size());

    return t_vec;
}

void Region::remove_trajectory(Trajectory* tr)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Region::removeTrajectory() called\n";
#endif

#ifdef DBG_ASSERT
    // when removing a trajectory, it should be already run out of segments
    if (!tr->segments.empty())
    {
        cout << "Big problem";
    }

    // when removing a trajectory, it should be already run out of tips (i.e. notification list empty)
    if (!tr->notificationList.empty())
    {
        cout << "oh dear\n";
    }
#endif

    // decrease total number of intersections in the geometry
    geometry->system->count_intersections -= 2 * (tr->intersections.size());

    // skip the first intersection iterator, as it is a garbage
    IntersectionItr is(++tr->intersections.begin());

    Trajectory* other_tr(nullptr);

    while (is != tr->intersections.end())
    {
        other_tr = is->second.other_trajectory;

        // invalidate all the intersetion for the other trajectory
        other_tr->invalidate_intersection(is->second.mirror);

        // erase all the intersections for the other trajectory
        other_tr->intersections.erase(is->second.mirror);

        is++;
    }

    // now clear all itersections of the trajectory
    tr->intersections.clear();

    // important to inform the linked (previous/next) trajectories about this removal
    if (tr->prev_tr.trajectory != nullptr)
    {
        if (tr->prev_tr.dir == BACKWARD)
        {
            tr->prev_tr.trajectory->next_tr.trajectory = nullptr;
        }
        else
        {
            tr->prev_tr.trajectory->prev_tr.trajectory = nullptr;
        }

        tr->prev_tr.trajectory = nullptr;
    }

    if (tr->next_tr.trajectory != nullptr)
    {
        if (tr->next_tr.dir == ::FORWARD)
        {
            tr->next_tr.trajectory->prev_tr.trajectory = nullptr;
        }
        else
        {
            tr->next_tr.trajectory->next_tr.trajectory = nullptr;
        }
        tr->next_tr.trajectory = nullptr;
    }

    // now remove the trajectory
    trajectories.remove(tr);
}

double Region::optical_length()
{
    Trajectory* tr(trajectories.first());
    double sum(0.0);

    // calculate total optical length present in this region
    while (tr != nullptr)
    {
        sum += tr->covered_length();
        tr = tr->next();
    }
    return sum;
}

RegionMTTipTag Region::register_on_region(MTTip* p_tip, TipType tiptype, MTType mttype)
{
    // update length
    update_region_length();

    // minus tip
    if (tiptype == T_MINUS)
    {
        return minus_tip_list.insert(minus_tip_list.end(), p_tip);
    }

    // plus tip
    else
    {
        // insert growing MT
        if (mttype == MT_GROWING)
        {
            geometry->system->growing_tips_reg[face_tag]++;
            return growing_plus_tip_list.insert(growing_plus_tip_list.end(), p_tip);
        }

        // insert shrinking MT
        else
        {
            return shrinking_plus_tip_list.insert(shrinking_plus_tip_list.end(), p_tip);
        }
    }
}

void Region::unregister_from_region(RegionMTTipTag tag, TipType tiptype, MTType mttype)
{
    // update length
    update_region_length();

    // minus tip
    if (tiptype == T_MINUS)
    {
        minus_tip_list.erase(tag);
    }

    // plus tip
    else
    {
        // erase growing MT
        if (mttype == MT_GROWING)
        {
            geometry->system->growing_tips_reg[face_tag]--;
            growing_plus_tip_list.erase(tag);
        }

        // erase shrinking MT
        else
        {
            shrinking_plus_tip_list.erase(tag);
        }
    }
    return;
}

void Region::update_region_length(bool force_update)
{
    // can go wrong on position cache refresh, a cleaner solution required
    if ((previous_update_tag == geometry->system->current_time_tag) && (!force_update))
    {
        return;
    }

    total_length += geometry->system->time_queue.progression(previous_update_tag)
                    * (geometry->system->p.v_min * (shrinking_plus_tip_list.size())
                       - geometry->system->p.v_tm * (minus_tip_list.size()))
                    + geometry->system->v_plus_queue.progression(previous_update_tag)
                      * (geometry->system->p.v_plus * (growing_plus_tip_list.size()));

    previous_update_tag = geometry->system->current_time_tag;

    return;
}
