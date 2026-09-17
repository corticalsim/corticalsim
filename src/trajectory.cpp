#include "trajectory.hpp"
#include "region.hpp"
#include "microtubule.hpp"
#include "segment.hpp"
#include "system.hpp"

Trajectory::Trajectory(SurfaceVector baseVec, vector<PointATedge> endP, double l):
    base(baseVec),
    length(l),
    prevTr(0, ::forward, NULL),
    nextTr(0, ::forward, NULL),
    prevTrCosAngle(1.0),
    nextTrCosAngle(1.0)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory created\n";
#endif

    // assign end points of the trajectory
    for (int i = 0; i < 2; i++)
    {
        endPoint.push_back(endP[i]);
    }

    // increase total number of trajectories in the geometry
    base.region->geometry->system->countTrajectories++;

    // initialize the intersection list of the trajectory
    intersections.insert(pair<double, Intersection>(-1.0, Intersection()));

    return;
}

Trajectory::~Trajectory()
{
    // decrease total number of trajectories in the geometry
    base.region->geometry->system->countTrajectories--;
    return;
}

bool Trajectory::integrityCheck()
{
    bool valid(true);

    if (prevTr.trajectory != NULL)
    {
        if (prevTr.dir == backward)
        {
            if (prevTr.trajectory->nextTr.trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
        else
        {
            if (prevTr.trajectory->prevTr.trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
    }
    if (nextTr.trajectory != NULL)
    {
        if (nextTr.dir == ::forward)
        {
            if (nextTr.trajectory->prevTr.trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
        else
        {
            if (nextTr.trajectory->nextTr.trajectory != this)
            {
                valid = false;
                cout << "Integrity: trajectory dependency chain broken.\n";
            }
        }
    }

    IntersectionItr is(wallBegin());

    for (size_t isCount = 1; isCount < intersections.size(); isCount++)
    {
        is++;
        if (is->second.mirror->second.mirror != is)
        {
            cout << "Integrity: intersection mirroring broken.\n";
            valid = false;
        }
    }

    list<Segment*>::iterator segItr(segments.begin());

    while (segItr != segments.end())
    {
        is = (**segItr).startItr;
        if ((**segItr).dir == ::forward)
        {
            while (++is != (**segItr).endItr)
            {
                is->second.occupancy--;
            }
        }
        else
        {
            while (--is != (**segItr).endItr)
            {
                is->second.occupancy--;
            }
        }
        segItr++;
    }
    // now, check whether all occupancy counts are zero
    is = intersections.begin();
    while (++is != intersections.end())
    {
        if (is->second.occupancy != 0)
        {
            cerr << "Detected occupancy " << is->second.occupancy << " (should be zero).\n";
            valid = false;
        }
    }

    segItr = segments.begin();
    while (segItr != segments.end())
    {
        is = (**segItr).startItr;
        if ((**segItr).dir == ::forward)
        {
            while (++is != (**segItr).endItr)
            {
                is->second.occupancy++;
            }
        }
        else
        {
            while (--is != (**segItr).endItr)
            {
                is->second.occupancy++;
            }
        }
        segItr++;
    }

    if (!valid)
    {
        cerr << "Invalid occupancy count found during integrity check.\n";
    }

    return valid;
}

double Trajectory::segmentLength()
{
    double sumLength(0.0);
    list<Segment*>::iterator seg(segments.begin());

    // calculate total length of all the segments present in the trajectory
    while (seg != segments.end())
    {
        sumLength += (**seg).length();
        seg++;
    }
    return sumLength;
}

double Trajectory::coveredLength()
{
    list<Segment*>::iterator seg(segments.begin());
    list<double> startPoints;
    list<double> endPoints;

    // store the start and end point locations of segments in two seperate list
    while (seg != segments.end())
    {
        startPoints.push_back(min((**seg).start, (**seg).end));
        endPoints.push_back(max((**seg).start, (**seg).end));
        seg++;
    }

    // make sure the last value of start-point-list is maximum of all values from both lists
    startPoints.push_back(VERY_LARGE);

    // sort the distances
    startPoints.sort();
    endPoints.sort();

    double previousPos(0.0);
    double currentPos(0.0);
    int occupancy(0);
    double cLength(0.0);
    list<double>::iterator stepup(startPoints.begin());
    list<double>::iterator stepdown(endPoints.begin());

    while (stepdown != endPoints.end())
    {
        currentPos = min(*stepup, *stepdown);

        // sum up only once and if the trajectory is covered, i.e. avoid multiple addition for bundles
        if (occupancy >= 1)
        {
            cLength += currentPos - previousPos;
        }

        previousPos = currentPos;

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

    return cLength;
}

void Trajectory::invalidateIntersection(IntersectionItr& oldIs)
{
    list<MTTip*>::iterator tip(notificationList.begin());

    // notify all tips on the current trajectory that a previous insertion is invalidated
    while (tip != notificationList.end())
    {
        (**tip).notifyRemove(oldIs);
        tip++;
    }

    return;
}

void Trajectory::newIntersection(IntersectionItr& newIs)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory::newIntersection() called.\n";
#endif

    list<Segment*>::iterator seg(segments.begin());

    int occupancy(0);

    // find whether the new intersection is a cross-intersection for the segment
    while (seg != segments.end())
    {
        (**seg).mt->updateLength();

        if ((**seg).crossesIntersection(newIs))
        {
            occupancy++;
        }

        seg++;
    }

    // assign complemetary occupancy for the new intersection
    newIs->second.occupancy = occupancy;

    // notify all tips on the current trajectory that a new insertion has happened
    list<MTTip*>::iterator tip(notificationList.begin());
    while (tip != notificationList.end())
    {
        (**tip).notifyInsert(newIs);
        tip++;
    }

    return;
}

int Trajectory::differenceSign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2)
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
        if (itr1 == wallEnd())
        {
            return 1;
        }

        do
        {
            itr2++;
        } while ((itr2 != wallEnd()) && (itr1 != itr2));

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

TrjSegmentTag Trajectory::insertSegment(Segment* s)
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
    return segments.insert(segments.end(), s);
}

void Trajectory::removeSegment(TrjSegmentTag s)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory::removeSegment() called.\n";
#endif

    // remove segment
    segments.erase(s);

    // if the associated trajectory has no segement and no tips, then remove the whole trajectory
    conditionalRemove();

    return;
}

TrjMTTipTag Trajectory::registerForNotifications(MTTip* pTip)
{
    // a new tip arrived on the trajectory, add it to the notification list
    return notificationList.insert(notificationList.end(), pTip);
}

void Trajectory::unregisterForNotifications(TrjMTTipTag tag)
{
    // a tip is leaving the trajectory, remove it from the notification list
    notificationList.erase(tag);

    // if the associated trajectory has no tips and no segments, then remove the whole trajectory
    conditionalRemove();
    return;
}

TrajectoryVector Trajectory::nextTrajectory(Direction dir)
{
    if (dir == ::forward)
    {
        // direction forward and a next trajectory is waiting to take the tip on board
        if (nextTr.trajectory != NULL)
        {
#ifdef DBG_GEOMETRY
            cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from "
                 << RegionTypeText[base.region->type] << " to " << RegionTypeText[nextTr.trajectory->base.region->type]
                 << "]\n";
#endif

            return nextTr;
        }
    }

    else
    {
        // direction backward and a previous trajectory is waiting to take the on board
        if (prevTr.trajectory != NULL)
        {
#ifdef DBG_GEOMETRY
            cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from "
                 << RegionTypeText[base.region->type] << " to " << RegionTypeText[prevTr.trajectory->base.region->type]
                 << "]\n";
#endif

            return prevTr;
        }
    }

    // there is no next/previous trajectory to take the tip, so extend (create a new) trajectory
    return base.region->geometry->extendTrajectory(this, dir);
}

void Trajectory::conditionalRemove()
{
    // if there is no segments and no tips, associted with a trajectory, then remove it rigth way
    if (segments.empty() && (notificationList.size() == 0))
    {
        base.region->removeTrajectory(this);
    }

    return;
}
