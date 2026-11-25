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
    t_vec.trajectory = trajectories.create(new_base, end_point, t_length);

    // make a list of intersections of the created trajectory with all other trajectories
    make_intersection_list(t_vec.trajectory);

    // increase total number of intersections accordingly
    geometry->system->countIntersections += 2 * (t_vec.trajectory->intersections.size());

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
    geometry->system->countIntersections -= 2 * (tr->intersections.size());

    // skip the first intersection iterator, as it is a garbage
    IntersectionItr is(++tr->intersections.begin());

    Trajectory* other_tr(nullptr);

    while (is != tr->intersections.end())
    {
        other_tr = is->second.otherTrajectory;

        // invalidate all the intersetion for the other trajectory
        other_tr->invalidate_intersection(is->second.mirror);

        // erase all the intersections for the other trajectory
        other_tr->intersections.erase(is->second.mirror);

        is++;
    }

    // now clear all itersections of the trajectory
    tr->intersections.clear();

    // important to inform the linked (previous/next) trajectories about this removal
    if (tr->prevTr.trajectory != nullptr)
    {
        if (tr->prevTr.dir == BACKWARD)
        {
            tr->prevTr.trajectory->nextTr.trajectory = nullptr;
        }
        else
        {
            tr->prevTr.trajectory->prevTr.trajectory = nullptr;
        }

        tr->prevTr.trajectory = nullptr;
    }

    if (tr->nextTr.trajectory != nullptr)
    {
        if (tr->nextTr.dir == ::FORWARD)
        {
            tr->nextTr.trajectory->prevTr.trajectory = nullptr;
        }
        else
        {
            tr->nextTr.trajectory->nextTr.trajectory = nullptr;
        }
        tr->nextTr.trajectory = nullptr;
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
        return minusTipList.insert(minusTipList.end(), p_tip);
    }

    // plus tip
    else
    {
        // insert growing MT
        if (mttype == MT_GROWING)
        {
            geometry->system->growingTipsReg[faceTag]++;
            return growingPlusTipList.insert(growingPlusTipList.end(), p_tip);
        }

        // insert shrinking MT
        else
        {
            return shrinkingPlusTipList.insert(shrinkingPlusTipList.end(), p_tip);
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
        minusTipList.erase(tag);
    }

    // plus tip
    else
    {
        // erase growing MT
        if (mttype == MT_GROWING)
        {
            geometry->system->growingTipsReg[faceTag]--;
            growingPlusTipList.erase(tag);
        }

        // erase shrinking MT
        else
        {
            shrinkingPlusTipList.erase(tag);
        }
    }
    return;
}

void Region::update_region_length(bool force_update)
{
    // can go wrong on position cache refresh, a cleaner solution required
    if ((previousUpdateTag == geometry->system->currentTimeTag) && (!force_update))
    {
        return;
    }

    totalLength
    += geometry->system->timeQueue.progression(previousUpdateTag)
       * (geometry->system->p.vMin * (shrinkingPlusTipList.size()) - geometry->system->p.vTM * (minusTipList.size()))
       + geometry->system->vPlusQueue.progression(previousUpdateTag)
         * (geometry->system->p.vPlus * (growingPlusTipList.size()));

    previousUpdateTag = geometry->system->currentTimeTag;

    return;
}
