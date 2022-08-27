#include "corticalSimReal.h"

void Geometry::callTranslator(SurfaceVector& sVec,int regionIndex,int regionIndexNeigh)
{
    int Index1(regions[regionIndex]->sideRevMap[regionIndexNeigh]);

    // use affine transformation of triangle to map a point from one region to other region 
    Vector2d V1(sVec.x,sVec.y),V2;
    V2 = regions[regionIndex]->side[Index1].A*V1 + regions[regionIndex]->side[Index1].b;
    sVec.x = V2(0,0);
    sVec.y = V2(1,0);

    sVec.angle+=regions[regionIndex]->side[Index1].xyRotationAngle;

    return;
}

TrajectoryVector Geometry::createTrajectory(const SurfaceVector& sVec)
{
    // create a new trajectory with give trajectory vector
    return sVec.region->insertTrajectory(sVec);
}

TrajectoryVector Geometry::createAndLinkTrajectory(const SurfaceVector& newBase, Trajectory* oldtr, Direction olddir, double cosAngle,double pCat)
{
    // first create a new trajectory vector (trajectory)
    TrajectoryVector tVec = createTrajectory(newBase);

    // declare a dummy trajectory vector 
    TrajectoryVector thisVec;

    // copy old trajectory vector
    thisVec.trajectory = oldtr;

    // old trajectory vector has direction forward, set the new trajectory vector as next 
    if (olddir == ::forward)
    {
        oldtr->nextTr = tVec;
        oldtr->nextTrCosAngle = cosAngle;
	oldtr->nextTrpCat = pCat;

        // get the direction and position of the dummy trajectory vector 
        thisVec.dir = backward;
        thisVec.pos = oldtr->length;
    }

    // old trajectory vector has direction backward, set the new trajectory vector as previous
    else
    {
        oldtr->prevTr = tVec;
        oldtr->prevTrCosAngle = cosAngle;
	oldtr->prevTrpCat = pCat;
        
        // get the direction and position of the dummy trajectory vector 
        thisVec.dir = ::forward;
        thisVec.pos = 0;
    }

    // new trajectory vector has direction forward, set the dummy trajectory vector as previous 
    if (tVec.dir == ::forward)
    {
        tVec.trajectory->prevTr = thisVec;
        tVec.trajectory->prevTrCosAngle = cosAngle;
        tVec.trajectory->prevTrpCat = pCat;
    }

    // new trajectory vector has direction backward, set the dummy trajectory vector as next
    else
    {
        tVec.trajectory->nextTr = thisVec;
        tVec.trajectory->nextTrCosAngle = cosAngle;
	tVec.trajectory->nextTrpCat = pCat;
    }

    return tVec;
}

SurfaceVector Geometry::randomSurfaceVector()
{
	int regionSelector = 0;
	double rArea = 0;
        int domainType;
                       
        if (regions.size() > 1)
	{
		// nucleation on PPB
		if(system->p.PPBkNucFraction > 0)
		{
			double randGen = system->randomGen.rand();
			domainType = (randGen < system->p.PPBkNucFraction) ? 1 : 2;
		}

		// nucleation outside PPB
		else
		domainType = 0;
                
		// select a region 
		rArea = system->randomGen.randDblExc(patchArea[domainType]);

		if (rArea < 0.5*patchArea[domainType])
		{
			regionSelector = 0;
			while (rArea > regions[RegionsIndex[domainType][regionSelector]]->area)
			{
				rArea -= regions[RegionsIndex[domainType][regionSelector++]]->area;

                		#ifdef DBG_ACID_TEST
				if (regionSelector == regions.size())
				{
					cerr << "DBG/ASSERT: ERROR: Ran out of surfaces for random surface vector creation.\n";
					exit(-1);
				}
                		#endif
			}
		}

		else
		{
			regionSelector = RegionsIndex[domainType].size()-1;
			rArea = patchArea[domainType] - rArea;
			while (rArea > regions[RegionsIndex[domainType][regionSelector]]->area)
			{
				rArea -= regions[RegionsIndex[domainType][regionSelector--]]->area;

                		#ifdef DBG_ACID_TEST
				if (regionSelector == -1)
				{
					cerr << "DBG/ASSERT: ERROR: Ran out of surfaces for random surface vector creation.\n";
					exit(-1);
				}
                		#endif
			}
		}
	}

	return regions[RegionsIndex[domainType][regionSelector]]->randomSurfaceVector();
}

bool Geometry::integrityCheck()
{
    int ridx(0);
    bool valid(true);
    int plusGrowCount(0);
    int plusShrinkCount(0);
    int minusCount(0);
    Trajectory* tptr(NULL);

    for (ridx = 0; ridx < regions.size(); ridx++)
    {
        plusGrowCount += regions[ridx]->growingPlusTipList.size();
        plusShrinkCount += regions[ridx]->shrinkingPlusTipList.size();
        minusCount += regions[ridx]->minusTipList.size();
        tptr = regions[ridx]->trajectories.first();
        while (tptr != NULL)
        {
            if (!tptr->integrityCheck())
                valid = false;

            tptr = tptr->next();
        }
    }

    if (minusCount != plusGrowCount + plusShrinkCount)
    {
        cerr << "Integrity check failed: unequal number of plus and minus ends registered on regions\n";
        valid = false;
    }

    if (plusGrowCount != system->growing_mts.size())
    {
        cerr << "Integrity check failed: number of registered growing tips does not equal the number of growing microtubules.\n";
        valid = false;
    }

    if (plusShrinkCount != system->shrinking_mts.size())
    {
        cerr << "Integrity check failed: number of registered shrinking tips does not equal the number of shrinking microtubules.\n";
        valid = false;
    }

    int totalPlusGrowCount(0);
    for(int dTag = 0; dTag <= system->p.faceNumber; dTag++)
    {
        totalPlusGrowCount+= system->growingTipsReg[dTag];
    }

    if (plusGrowCount != totalPlusGrowCount)
    {
        cerr << "Integrity check failed: number of growing tips does not equal cached values\n";
        valid = false;
    }

    return valid;
}


double Geometry::opticalLength()
{
    int ridx(0);
    double sum(0.0);

    // calculate total optical length present in the geometry
    for (ridx = 0; ridx < regions.size(); ridx++)
    {
        sum += regions[ridx]->opticalLength();
    }
    return sum;
}

int Geometry::trajectoryCount()
{
    int ridx(0);
    int sum(0);

    // calculate total number of trajectory present in the geometry
    for (ridx = 0; ridx < regions.size(); ridx++)
    {
        sum += regions[ridx]->trajectories.size();
    }

    return sum;
}

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
    getTrajectoryCoordinates(newBase,tLength,endPoint,tVec);

    // create the trajectory
    tVec.trajectory = trajectories.create(newBase,endPoint,tLength);

    // make a list of intersections of the created trajectory with all other trajectories 
    makeIntersectionList(tVec.trajectory);

    // increase total number of intersections accordingly 
    geometry->system->countIntersections += 2*(tVec.trajectory->intersections.size());

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
    geometry->system->countIntersections -= 2* (tr->intersections.size());

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
            tr->prevTr.trajectory->nextTr.trajectory = NULL;
        else
            tr->prevTr.trajectory->prevTr.trajectory = NULL;

        tr->prevTr.trajectory = NULL;
    }

    if (tr->nextTr.trajectory != NULL)
    {
        if (tr->nextTr.dir == ::forward)
            tr->nextTr.trajectory->prevTr.trajectory = NULL;
        else
            tr->nextTr.trajectory->nextTr.trajectory = NULL;
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
        return minusTipList.insert(minusTipList.end(),pTip);

    // plus tip
    else
    {
        // insert growing MT
        if (mttype == mt_growing)
        {
	    geometry->system->growingTipsReg[faceTag]++;
            return growingPlusTipList.insert(growingPlusTipList.end(),pTip);
        }

        // insert shrinking MT 
        else
            return shrinkingPlusTipList.insert(shrinkingPlusTipList.end(),pTip);
    }
}

void Region::unregisterFromRegion(RegionMTTipTag tag, TipType tiptype, MTType mttype)
{
    // update length
    updateRegionLength();

    // minus tip
    if (tiptype == t_minus)
        minusTipList.erase(tag);

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
            shrinkingPlusTipList.erase(tag);
    }
    return;
}

Trajectory::Trajectory(SurfaceVector baseVec,vector<PointATedge> endP,double l) :
    base(baseVec),
    length(l),
    prevTr(0,::forward,NULL),
    nextTr(0,::forward,NULL),
    prevTrCosAngle(1.0),
    nextTrCosAngle(1.0)
{
    #ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Trajectory created\n";
    #endif

    // assign end points of the trajectory
    for(int i=0;i<2;i++)
    endPoint.push_back(endP[i]);

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
    int isCount(0);

    for (isCount = 1; isCount < intersections.size(); isCount++)
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
                is->second.occupancy--;
        }
        else
        {
            while (--is != (**segItr).endItr)
                is->second.occupancy--;
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
                is->second.occupancy++;
        }
        else
        {
            while (--is != (**segItr).endItr)
                is->second.occupancy++;
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
            cLength += currentPos - previousPos;

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

#ifndef NO_INLINE
inline
#endif

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
    // returns the sign of 'pos1 - pos2'. If the two positions are equal, the iterators are used instead to determine order.
    if(pos1 - pos2 > ZERO_CUTOFF)
        return 1;

    else if(pos1 - pos2 < -ZERO_CUTOFF)
        return -1;

    else
    {
        if(itr1 == itr2)
            return 0;
        if(itr1 == wallEnd())
            return 1;

        do
        {
            itr2++;
        }
        while ((itr2 != wallEnd()) && (itr1 != itr2));

        if(itr1 == itr2)
            return 1;
        else
            return -1;
    }
}

TrjSegmentTag Trajectory::insertSegment(Segment* s)
{
    #ifdef DBG_ASSERT
    // check whether the size of the inserted segment equals zero, for insertion of non-zero segments, take care to update the occupancy numbers
    if (abs(s->end - s->start) > ZERO_CUTOFF)
    {
        cerr << "DBG/ASSERT: ERROR: Not permitted to insert segment with non-zero length.\n";
        exit(-1);
    }
    #endif

    // insert a new segment to the trajectory 
    return segments.insert(segments.end(),s);
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
    return notificationList.insert(notificationList.end(),pTip);
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
            cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from " << RegionTypeText[base.region->type] << \
                 " to " << RegionTypeText[nextTr.trajectory->base.region->type] << "]\n";
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
            cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from " << RegionTypeText[base.region->type] << \
                 " to " << RegionTypeText[prevTr.trajectory->base.region->type] << "]\n";
            #endif

            return prevTr;
        }
    }

    // there is no next/previous trajectory to take the tip, so extend (create a new) trajectory 
    return base.region->geometry->extendTrajectory(this,dir);
}
