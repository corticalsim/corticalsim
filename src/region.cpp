#include "region.hpp"
#include "geometry.hpp"
#include "system.hpp"

TrajectoryVector Region::insertTrajectory(const SurfaceVector& sVec)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Region::insertTrajectory() called\n";
#endif

    TrajectoryVector tVec;
    double tLength(0.);
    SurfaceVector newBase(sVec);

    vector<PointATedge> endPoint;

    // gather coordinates of the trajectory to be created
    getTrajectoryCoordinates(newBase, tLength, endPoint, tVec);

    // create the trajectory
    tVec.trajectory = trajectories.create(newBase, endPoint, tLength);

    // make a list of intersections of the created trajectory with all other trajectories
    makeIntersectionList(tVec.trajectory);

    // increase total number of intersections accordingly
    geometry->system->countIntersections += 2 * (tVec.trajectory->intersections.size());

    return tVec;
}

void Region::removeTrajectory(Trajectory* tr)
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

    Trajectory* otherTr(NULL);

    while (is != tr->intersections.end())
    {
        otherTr = is->second.otherTrajectory;

        // invalidate all the intersetion for the other trajectory
        otherTr->invalidateIntersection(is->second.mirror);

        // erase all the intersections for the other trajectory
        otherTr->intersections.erase(is->second.mirror);

        is++;
    }

    // now clear all itersections of the trajectory
    tr->intersections.clear();

    // important to inform the linked (previous/next) trajectories about this removal
    if (tr->prevTr.trajectory != NULL)
    {
        if (tr->prevTr.dir == backward)
        {
            tr->prevTr.trajectory->nextTr.trajectory = NULL;
        }
        else
        {
            tr->prevTr.trajectory->prevTr.trajectory = NULL;
        }

        tr->prevTr.trajectory = NULL;
    }

    if (tr->nextTr.trajectory != NULL)
    {
        if (tr->nextTr.dir == ::forward)
        {
            tr->nextTr.trajectory->prevTr.trajectory = NULL;
        }
        else
        {
            tr->nextTr.trajectory->nextTr.trajectory = NULL;
        }
        tr->nextTr.trajectory = NULL;
    }

    // now remove the trajectory
    trajectories.remove(tr);
}

double Region::opticalLength()
{
    Trajectory* tr(trajectories.first());
    double sum(0.0);

    // calculate total optical length present in this region
    while (tr != NULL)
    {
        sum += tr->coveredLength();
        tr = tr->next();
    }
    return sum;
}

RegionMTTipTag Region::registerOnRegion(MTTip* pTip, TipType tiptype, MTType mttype)
{
    // update length
    updateRegionLength();

    // minus tip
    if (tiptype == t_minus)
    {
        return minusTipList.insert(minusTipList.end(), pTip);
    }

    // plus tip
    else
    {
        // insert growing MT
        if (mttype == mt_growing)
        {
            geometry->system->growingTipsReg[faceTag]++;
            return growingPlusTipList.insert(growingPlusTipList.end(), pTip);
        }

        // insert shrinking MT
        else
        {
            return shrinkingPlusTipList.insert(shrinkingPlusTipList.end(), pTip);
        }
    }
}

void Region::unregisterFromRegion(RegionMTTipTag tag, TipType tiptype, MTType mttype)
{
    // update length
    updateRegionLength();

    // minus tip
    if (tiptype == t_minus)
    {
        minusTipList.erase(tag);
    }

    // plus tip
    else
    {
        // erase growing MT
        if (mttype == mt_growing)
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

void Region::updateRegionLength(bool forceUpdate)
{
    // can go wrong on position cache refresh, a cleaner solution required
    if ((previousUpdateTag == geometry->system->currentTimeTag) && (!forceUpdate))
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
