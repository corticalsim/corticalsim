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

    vector<PointATedge> end_point;

    // gather coordinates of the trajectory to be created
    get_trajectory_coordinates(new_base, t_length, end_point, t_vec);

    // create the trajectory
    t_vec.m_trajectory = m_trajectories.create(new_base, end_point, t_length);

    // make a list of intersections of the created trajectory with all other trajectories
    make_intersection_list(t_vec.m_trajectory);

    // increase total number of intersections accordingly
    m_geometry->m_system->m_count_intersections += 2 * (t_vec.m_trajectory->m_intersections.size());

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
    m_geometry->m_system->m_count_intersections -= 2 * (tr->m_intersections.size());

    // skip the first intersection iterator, as it is a garbage
    IntersectionItr is(++tr->m_intersections.begin());

    Trajectory* other_tr(nullptr);

    while (is != tr->m_intersections.end())
    {
        other_tr = is->second.m_other_trajectory;

        // invalidate all the intersetion for the other trajectory
        other_tr->invalidate_intersection(is->second.m_mirror);

        // erase all the intersections for the other trajectory
        other_tr->m_intersections.erase(is->second.m_mirror);

        is++;
    }

    // now clear all itersections of the trajectory
    tr->m_intersections.clear();

    // important to inform the linked (previous/next) trajectories about this removal
    if (tr->m_prev_tr.m_trajectory != nullptr)
    {
        if (tr->m_prev_tr.m_dir == BACKWARD)
        {
            tr->m_prev_tr.m_trajectory->m_next_tr.m_trajectory = nullptr;
        }
        else
        {
            tr->m_prev_tr.m_trajectory->m_prev_tr.m_trajectory = nullptr;
        }

        tr->m_prev_tr.m_trajectory = nullptr;
    }

    if (tr->m_next_tr.m_trajectory != nullptr)
    {
        if (tr->m_next_tr.m_dir == ::FORWARD)
        {
            tr->m_next_tr.m_trajectory->m_prev_tr.m_trajectory = nullptr;
        }
        else
        {
            tr->m_next_tr.m_trajectory->m_next_tr.m_trajectory = nullptr;
        }
        tr->m_next_tr.m_trajectory = nullptr;
    }

    // now remove the trajectory
    m_trajectories.remove(tr);
}

double Region::optical_length()
{
    Trajectory* tr(m_trajectories.first());
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
        return m_minus_tip_list.insert(m_minus_tip_list.end(), p_tip);
    }

    // plus tip
    else
    {
        // insert growing MT
        if (mttype == MT_GROWING)
        {
            m_geometry->m_system->m_growing_tips_reg[m_face_tag]++;
            return m_growing_plus_tip_list.insert(m_growing_plus_tip_list.end(), p_tip);
        }

        // insert shrinking MT
        else
        {
            return m_shrinking_plus_tip_list.insert(m_shrinking_plus_tip_list.end(), p_tip);
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
        m_minus_tip_list.erase(tag);
    }

    // plus tip
    else
    {
        // erase growing MT
        if (mttype == MT_GROWING)
        {
            m_geometry->m_system->m_growing_tips_reg[m_face_tag]--;
            m_growing_plus_tip_list.erase(tag);
        }

        // erase shrinking MT
        else
        {
            m_shrinking_plus_tip_list.erase(tag);
        }
    }
    return;
}

void Region::update_region_length(bool force_update)
{
    // can go wrong on position cache refresh, a cleaner solution required
    if ((m_previous_update_tag == m_geometry->m_system->m_current_time_tag) && (!force_update))
    {
        return;
    }

    m_total_length
    += m_geometry->m_system->m_time_queue.progression(m_previous_update_tag)
       * (m_geometry->m_system->m_p.m_v_min * (m_shrinking_plus_tip_list.size()) - m_geometry->m_system->m_p.m_v_tm * (m_minus_tip_list.size()))
       + m_geometry->m_system->m_v_plus_queue.progression(m_previous_update_tag)
         * (m_geometry->m_system->m_p.m_v_plus * (m_growing_plus_tip_list.size()));

    m_previous_update_tag = m_geometry->m_system->m_current_time_tag;

    return;
}
