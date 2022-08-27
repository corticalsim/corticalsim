#include "corticalSimReal.h"

Microtubule::Microtubule(System* s, TrajectoryVector tv, bool initialize)
    : system(s),
      plus(this, tv, s->p.vPlus, &(s->vPlusQueue), s->p.vPlus),

      // reverse direction
      minus(this, tv.flipped(), -s->p.vTM, &(s->timeQueue), -s->p.vTM),	
      disappearEvent(this, &(s->timeQueue), - s->p.vMin + s->p.vTM),
      type(mt_growing),
      nucleationTime(s->systemTime + s->systemTimeOffset),
      previousUpdateTag(s->currentTimeTag)
{
    #ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule created.\n";
    #endif

    if (initialize)
    {
        // create a single segment MT
        segments.create(this,tv);

        // intialize the plus tip of MT
        plus.initialize();

       // initialize the minus tip of MT
        minus.initialize();

        // microtubule is disappearing (segment length finishing)
        setDisappearEvent();
    }
    return;
}

Microtubule::~Microtubule()
{
    #ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule destroyed. address=" << this << "\n";
    #endif

    #ifdef DBG_ASSERT
    if ((segments.first()->start - segments.last()->end) > ZERO_CUTOFF)
    {
        cerr << "DBG/ASSERT: Microtubule::~Microtubule: endpoints do not coincide/\n";
        exit(-1);
    }
    #endif

    // un-linking the MT from the trajectory(ies)
    plus.unlinkFromTrajectory();
    minus.unlinkFromTrajectory();

    return;
}

void Microtubule::setDisappearEvent()
{
    // push to disappear event: (plus.velocity - minus.velocity < 0)
    if ((segments.size() == 1) && (type == mt_shrinking)) 
        disappearEvent.pushOnQueue(segments.first()->length(), ev_disappear);

    // flush disappear event queue
    else
        disappearEvent.clear();
    return;
}

void Microtubule::handleEvent(const EventDescriptor* ei)
{
    // better update length before processing next event 
    updateLength();

    #ifdef DBG_ACID_TEST
    if (!integrityCheck())
    {
        cerr << "failed before event\n";
        exit(-1);
    }
    #endif

    switch(ei->type)
    {
        // wall collision event encountered 
    	case ev_wall:
        	wall();
        break;

	// MT-MT collision event encountered 
    	case ev_collision:
        	collision();
        break;

        // back track event encountered
    	case ev_backtrack:
        	if (ei == &(plus.event))
            	backtrack(&plus);
        	else
            	backtrack(&minus);
        break;

	// end of segment event encountered
    	case ev_end_of_segment:
        	if (ei == &(plus.event))
            	endOfSegment(&plus);
        	else
            	endOfSegment(&minus);
        break;

        // disappear event encountered
    	case ev_disappear:
        	harakiri();
        break;

    };

    #ifdef DBG_ACID_TEST
    if ((ei->type != ev_disappear) && (!integrityCheck()))
    {
        cerr << "failed after event\n";
        exit(-1);
    }
    #endif

    return;
}

void Microtubule::catastrophe()
{
    #ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule::catastrophe() called.\n";
    #endif

    #ifdef DBG_EXTRA_CHECK
    if (type != mt_growing)
    {
        cout << "ERROR: catastrophe called for non-growing MT\n";
        exit(-1);
    }
    #endif

    // update MT length
    updateLength();

    // growing MT will become shrinking MT, hence unregister the tip from growing-plus-tip list of this region 
    plus.trajectory->base.region->unregisterFromRegion(plus.regionTag,t_plus,mt_growing);

    // now, this is shrinking MT 
    type = mt_shrinking;

    // growing MT will become shrinking MT, hence register the tip on shrinking-plus-tip list of this region 
    plus.regionTag = plus.trajectory->base.region->registerOnRegion(&plus,t_plus,mt_shrinking);

    // replace plus velocity by minus velocity  
    plus.velocity = system->p.vMin;

    // 
    plus.event.reinitialize(&(system->timeQueue),system->p.vMin);

    // import the growing MT as shrinking MT 
    system->shrinking_mts.import(system->growing_mts, this);

    if (system->p.vMin < 0)
    {
        // velocity has changed sign, select a new event
        plus.advanceIntersection();
    }

    // schedule next event, (if this is triggered by a collision, and vMin > 0, increase occupancy)
    plus.determineEvent();

    // check possibility of a disappear event 
    setDisappearEvent();

    return;
}

void Microtubule::rescue()
{
    #ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule::rescue() called.\n";
    #endif

    #ifdef DBG_EXTRA_CHECK
    if (type != mt_shrinking)
    {
        cout << "ERROR: rescue called for non-shrinking MT\n";
        exit(-1);
    }
    #endif

    updateLength();

    // shrinking MT will become growing  MT, hence unregister the tip from shrinking-plus-tip list of this region 
    plus.trajectory->base.region->unregisterFromRegion(plus.regionTag,t_plus,mt_shrinking);

    // now, this is growing MT 
    type = mt_growing;

    // shrinking MT will become growing MT, hence register the tip on growing-plus-tip list of this region
    plus.regionTag = plus.trajectory->base.region->registerOnRegion(&plus,t_plus,mt_growing);

    // replace minue velocity by plus velocity 
    plus.velocity = system->p.vPlus;

    plus.event.reinitialize(&(system->vPlusQueue), system->p.vPlus);

    // import the shrinking MT as growing MT
    system->growing_mts.import(system->shrinking_mts, this);

    if (system->p.vMin < 0)
    {
        // velocity has changed sign, need to select a new event
        plus.advanceIntersection();
    }

    // schedule next event
    plus.determineEvent();

    // check possibility of a disappear event
    setDisappearEvent();

    return;
}

void Microtubule::sever(Segment* cutSeg, double cutPos)
{
    // sever MT at given position
    splitSegmentAtTrajPos(cutSeg->start + cutSeg->dir*cutPos, cutSeg);
    system->totalLengthSeveringCount++;
    return;
}

void Microtubule::severAtCross(IntersectionItr is, Segment* cutSeg)
{
    // WARNING: take care to cut the segment right *after* the intersection. Not doing so can lead to problems
    // when treadmilling is disabled. The microtubule end is then located *exactly* on the intersection and can be cut again and again,
    // leading to problems.

    // this time, cutPos and cutSeg are already known.
    double splitPos(is->first);
    // shift the cut ever so slightly away from the intersection, in the direction of the plus end.
    if (cutSeg->dir == ::forward)
        splitPos += ZERO_CUTOFF;
    else
        splitPos -= ZERO_CUTOFF;

    // only cut if the severing position is comfortably located within the segment interval
    if ((cutSeg->start - splitPos)*(cutSeg->end - splitPos) < ZERO_CUTOFF)
    {
        splitSegmentAtTrajPos(splitPos, cutSeg);
        system->totalIntersectionSeveringCount++;
    }
    return;
}

void Microtubule::splitSegmentAtTrajPos(double cutPos, Segment* cutSeg)
{
    // used in sever and severAtCross
    Microtubule* newMT(NULL);
    TrajectoryVector tv;

    tv.trajectory = cutSeg->trajectory;
    tv.dir = cutSeg->dir;
    tv.pos = cutPos;
    if (type == mt_growing)
        newMT = system->growing_mts.create(system, tv, false);
    else
        newMT = system->shrinking_mts.create(system, tv, false);

    newMT->nucleationTime = nucleationTime;
    newMT->type = type;

    // transfer segments
    newMT->segments.create(newMT,tv);
    newMT->segments.first()->end = cutSeg->end;
    newMT->segments.first()->endItr = cutSeg->endItr;
    newMT->segments.first()->nucleationTime = cutSeg->nucleationTime;
    newMT->segments.importSet(segments, cutSeg->next(), segments.last());
    Segment* tempSeg = newMT->segments.first()->next();
    while (tempSeg != NULL)
    {
        tempSeg->mt = newMT;
        tempSeg = tempSeg->next();
    }

    // now update the plus tip trajectory/position
    newMT->plus.trajectory = plus.trajectory;
    newMT->plus.velocity = plus.velocity;
    newMT->plus.event.reinitialize(plus.event.queue ,plus.velocity);
    newMT->plus.dir = plus.dir;
    newMT->plus.nextCollision = plus.nextCollision;
    newMT->plus.notificationTag = newMT->plus.trajectory->registerForNotifications(&(newMT->plus));
    newMT->plus.regionTag = newMT->plus.trajectory->base.region->registerOnRegion(&(newMT->plus),t_plus,type);

    // next, update the minus end of the new MT
    newMT->segments.first()->start = tv.pos;
    newMT->minus.notificationTag = newMT->minus.trajectory->registerForNotifications(&(newMT->minus));
    newMT->minus.regionTag = newMT->minus.trajectory->base.region->registerOnRegion(&(newMT->minus),t_minus,type);
    newMT->minus.locateIntersection();

    newMT->plus.determineEvent();
    newMT->minus.determineEvent();
    newMT->setDisappearEvent();

    // now, finish up the base microtubule
    segments.last()->end = cutPos;

    // always unregister from region (just to be sure)
    plus.trajectory->base.region->unregisterFromRegion(plus.regionTag,t_plus,type);
    if (type == mt_growing)
    {
        type = mt_shrinking;
        plus.velocity = system->p.vMin;
        plus.event.reinitialize(&(system->timeQueue),system->p.vMin);
        system->shrinking_mts.import(system->growing_mts, this);
    }

    if (plus.trajectory != tv.trajectory)
    {
        plus.trajectory->unregisterForNotifications(plus.notificationTag);
        plus.trajectory = tv.trajectory;
        plus.notificationTag = plus.trajectory->registerForNotifications(&(plus));
    }
    //...and register on region again
    plus.regionTag = plus.trajectory->base.region->registerOnRegion(&plus,t_plus,mt_shrinking);
    plus.dir = tv.dir;
    plus.nextCollision = newMT->minus.nextCollision;
    if (system->p.vMin * newMT->minus.velocity >= 0)
    {
        // both tips move in opposite directions [this should generally be the case]
        plus.advanceIntersection();
    }

    plus.determineEvent();
    if (segments.size() == 1)
    {
        minus.determineEvent();
    }
    setDisappearEvent();

    return;
}

void Microtubule::translatePositionMT2Segment(double& cutPos, Segment*& cutSeg)
{
    cutSeg = segments.first();
    while ((cutSeg != NULL) && (cutPos > cutSeg->length()))
    {
        cutPos -= cutSeg->length();
        cutSeg = cutSeg->next();
    }
    if (cutPos > cutSeg->length())
    {
        if (cutPos - cutSeg->length() < ZERO_CUTOFF)
        {
            cutPos = cutSeg->length() - ZERO_CUTOFF;
        }
        else
        {
            cerr << "ERROR: position beyond the end of the microtubule.\n";
            exit (-1);
        }
    }
    cutPos = cutSeg->start + cutSeg->dir*cutPos;
    return;
}

void Microtubule::wall()
{
    // only for plus end 

    TrajectoryVector tv;

    // set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
    segments.last()->end = plus.nextEventPos;
    tv = plus.trajectory->nextTrajectory(plus.dir);

    // next trajectory on 'forbidden zone', leads catastrophe
    if (tv.trajectory->base.region->faceTag == -1)
    {
        tv.trajectory->conditionalRemove();
        catastrophe();

        // no wall crossing
        return;
    }

    // edge catastrophe enabled 
    if (system->p.edgeCatastropheEnabled)
    {
        double cosAngle = (plus.dir == ::forward) ? plus.trajectory->nextTrCosAngle : plus.trajectory->prevTrCosAngle;

        // MT is moving non-parallel to the edege 
        if ((1.0 - cosAngle) > ZERO_CUTOFF)
        {          
            double pCat(0.0);
	    
	    // for sharp edged cube 
	    if(system->p.geometry == "cubeReal")
	    {
		if(abs(plus.trajectory->base.region->faceTag - tv.trajectory->base.region->faceTag) > 0)
                {
		    	if (abs(plus.trajectory->base.region->faceTag - tv.trajectory->base.region->faceTag) > 1)
			pCat = system->p.pCatSpecialEdgeMax;

		    	else
		    	pCat = system->p.pCatRegularEdgeMax;
		}
	    }

	    // other  
            else
	    pCat = (plus.dir == ::forward) ? plus.trajectory->nextTrpCat : plus.trajectory->prevTrpCat; 

            // when smooth catastrophe is enabled 
            if(system->p.edgeCatastropheSmooth)
            {
                pCat *= (1.0 - cosAngle);
            }

            // apply catastrophe 
            if (system->randomGen.rand() <= pCat)
            {
                tv.trajectory->conditionalRemove();
                catastrophe();

                // no wall crossing
                return;
            }
        }
    }

    // update boundarty crossing count 
    system->boundaryCrossingCount++;

    #ifdef DBG_ASSERT
    if 	(segments.last()->endItr != plus.nextCollision)
        cerr << "DBG/ASSERT: nextCollision does not match with endItr upon entering wall event.\n";
    #endif

    // wall crossing 
    segments.last()->endItr = plus.nextCollision;

    // create a new segment 
    segments.create(this, tv);
    segments.last()->startItr = tv.dir == ::forward ? tv.trajectory->wallBegin() : tv.trajectory->wallEnd();

    // switch trajectory 
    plus.switchTrajectory(tv.trajectory, tv.dir, tv.dir == ::forward ? tv.trajectory->wallBegin() : tv.trajectory->wallEnd());

    if (segments.size() == 2)
    {
        // MT with two segment of zero (=0) length, consider possibilty of a disappear event 
        minus.determineEvent();
        setDisappearEvent();
    }

    return;
}

void Microtubule::collision()
{
    // set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
    segments.last()->end = plus.nextEventPos;
    CollisionType type;
    Direction dir;

    // no MT(s) passing throught the next collision point, resulting a simple crossover 
    if (plus.nextCollision->second.mirror->second.occupancy == 0)
        type = ct_crossover;

    // MT(s) passing through the next collision point  
    else
    {
        double diffAngle = plus.trajectory->base.region->intersectionAngle(plus.trajectory, plus.nextCollision->second.otherTrajectory);

        if (diffAngle > 0.5*PI)
        {
            	dir = (plus.dir == ::forward) ? backward : ::forward;
            	diffAngle = PI - diffAngle;
        }
        else
        	dir = plus.dir;
	
	// find out type of collision event to be occured 
        type = system->collisionResult(diffAngle, plus.nextCollision->second.occupancy, plus.nextCollision->second.mirror->second.occupancy);
    }

    switch(type)
    {
	// zippering 
    	case ct_zipper:
        	system->totalZipperCount++;
        	zipper(dir);
        break;

	// crossover 
    	case ct_crossover:
        	system->totalCrossoverCount++;
        	crossover();
        break;

	// induced catastrophe
    	case ct_inducedCatastrophe:
        	system->totalInducedCatastropheCount++;
        	catastrophe();
        break;

    #ifdef DBG_EXTRA_CHECK
    default:
        cerr << "Unknown collision result type (" << type << ") used.\n";
        exit(-1);
    #endif

    };

    return;
}

void Microtubule::zipper(Direction dir)
{
    #ifdef DBG_ASSERT
    if (abs(plus.nextCollision->first - segments.last()->end) > ZERO_CUTOFF)
    {
        cout << "DBG/ASSERT: ERROR: end of segment does not match zippering location.[discrepancy: " << abs(plus.nextCollision->first - segments.last()->end) << "]\n";
        exit(-1);
    }
    #endif

    // get a trajectory vetor on the new trajectory to be associted after zippering
    TrajectoryVector tv(plus.nextCollision->second.mirror->first, dir, plus.nextCollision->second.otherTrajectory);

    // link the end of the last segment (containing the plus tip) to the next plus end collision
    segments.last()->endItr = plus.nextCollision;

    // create a new zero (=0) length segment on new position and direction (of course on new trajectory)
    segments.create(this,tv);

    // link the bottom of the last segment (containing the new minus tip) to the next plus end collision mirror
    segments.last()->startItr = plus.nextCollision->second.mirror;

    // now do the switching of trajectory from incoming to outgoing
    plus.switchTrajectory(plus.nextCollision->second.otherTrajectory, dir, plus.nextCollision->second.mirror);

    if (segments.size() == 2)
    {
        // number of segment is two (=2), HERE it means two zero length segment, prepare to schedule a disappear event 
        minus.determineEvent();
        setDisappearEvent();
    }

    return;
}

void Microtubule::endOfSegment(MTTip* tip)
{
    #ifdef DBG_ASSERT
	if (segments.size() == 1)
	{
		cout << "DBG/ASSERT: ERROR: end-of-segment event called for single-segmented microtubule.\n";
		//exit(-1);
	}
	#endif

	Direction newDir;
	Segment* newSeg(NULL);
	Segment* killSeg(NULL);

        // select plus tip associated segment 
	if (tip == &plus)
	{
		killSeg = segments.last();
		newSeg = segments.last()->previous();
		newDir = newSeg->dir;
	}

        // select minus tip associated segment 
	else
	{
		killSeg = segments.first();
		newSeg = segments.first()->next();
		if (newSeg->dir == ::forward)
			newDir = backward;
		else
			newDir = ::forward;
	}

    #ifdef DBG_ASSERT
    if (abs(killSeg->end - killSeg->start) > ZERO_CUTOFF)
    {
        cout << "DEBUG/ASSERT: ERROR: non-zero segment length on end of segment event\n";
        exit(-1);
    }
    #endif

    // kill the selected segment 
    segments.remove(killSeg);

    // tip switching to a new trajectory (in a different region) after encountering wall-begin/wall-end 
    if ((tip->nextCollision == tip->trajectory->wallBegin()) ||(tip->nextCollision == tip->trajectory->wallEnd()))
    {
        tip->switchTrajectory(newSeg->trajectory, newDir, newDir == ::forward ? newSeg->trajectory->wallEnd() : newSeg->trajectory->wallBegin());
        system->boundaryCrossingCount--;
    }

    // tip switching to a new trajectory (in the same region) 
    else
        tip->switchTrajectory(tip->nextCollision->second.otherTrajectory, newDir, tip->nextCollision->second.mirror);

    // MT length zero (=0) encountered
    if (segments.size()==1)
    {
        // current tip is a plus tip, so determine the minus tip associated event 
        if (tip == &(plus))
            minus.determineEvent();

        // current tip is a minus tip, so determine the plus tip associated event 
        else
            plus.determineEvent();

        // schedule a disappear event 
        setDisappearEvent();
    }

    return;
}

void Microtubule::backtrack(MTTip* tip)
{
    // asociated with depolymerizing plus/minus end 

    if (tip == &(plus))
        segments.last()->end = plus.nextEventPos;
    else
        segments.first()->start = minus.nextEventPos;

    // MT is back tracked from this site, resulting a decrease in occupancy 
    tip->nextCollision->second.occupancy--;

    #ifdef CROSS_SEV
    if ((tip->nextCollision->second.occupancy == 0) && (tip->nextCollision->second.mirror->second.occupancy > 0))
    {
        system->removeOccupiedIntersection(tip->nextCollision->second);
    }
    #endif

    #ifdef DBG_ASSERT
    if (tip->nextCollision->second.occupancy < 0)
    {
        cerr << "DBG/ASSERT: ERROR: negative intersection occupancy after ";
        if (tip == &(plus))
            cerr << "plus";
        else
            cerr << "minus";
        cerr << " end backtracking.\n";

        cout << "number of equal elements: " << tip->trajectory->intersections.count(tip->nextCollision->first) << "\n";;

        double temp = tip->nextCollision->first;

        do
        {
            cout << "pos=" << tip->mt->segments.last()->end << ", " << tip->mt->segments.first()->start << ", occupancy=" << tip->nextCollision->second.occupancy << "  ";
            tip->advanceIntersection();
            tip->determineEvent();
            cout << "step (size=" << tip->nextCollision->first - temp << ")\n";
        }
        while (tip->event.type != ev_end_of_segment);
        cout << "pos=" << tip->mt->segments.last()->end << ", " << tip->mt->segments.first()->start << "  \n";

        exit(-1);
    }
    #endif

    // schedule next event 
    tip->advanceIntersection();
    tip->determineEvent();

    return;
}

void Microtubule::crossover()
{
    #ifdef DBG_ASSERT
    if (plus.nextCollision->second.occupancy < 0)
    {
        cerr << "DBG/ASSERT: ERROR: negative intersection occupancy before crossover.\n";
        exit(-1);
    }
    #endif

    #ifdef CROSS_SEV
    if ((plus.nextCollision->second.occupancy == 0) && plus.nextCollision->second.mirror->second.occupancy > 0)
    {
        system->addOccupiedIntersection(plus.nextCollision);
    }
    #endif

    // skip any (active) event in the current intersection point, instead move the tip to the next intersection point 
    plus.nextCollision->second.occupancy++;
    plus.advanceIntersection();
    plus.determineEvent();

    return;
}

void Microtubule::harakiri()
{
    // MT disappear event encountered, remove the MT
    if (type == mt_growing)
        system->growing_mts.remove(this);
    else
        system->shrinking_mts.remove(this);

    return;
}

bool Microtubule::integrityCheck()
{
    bool valid(true);
    double maxDeviation(0.0);

    Segment* seg(NULL);
    IntersectionItr is;
    int count(0);

    // check connections between end-of segment iterators and match with positions (should be exact)
    if (segments.size() == 0)
    {
        // MT should contain at least one segment
        valid = false;
        cerr << "Microtubule contains no segments.\n";
    }

    else if (segments.size() != 1)
    {
        if (segments.size() == 0)
        {
            // MT should contain at least one segment
            cerr  << "microtubule without segments!\n";
            exit(-1);
        }

        seg = segments.first();
        if (seg == NULL)
        {
            // pointer to segment list can not be null
            cerr << "null pointer received\n";
            exit(-1);
        }

        // segment (some where) at the middle of a trajectory
        if ((seg->endItr != seg->trajectory->wallEnd()) && (seg->endItr != seg->trajectory->wallBegin()))
        {
            if (seg->endItr->second.mirror != seg->next()->startItr)
            {
                // next-start-iterator should match with the end-second-mirror-iterator 
                valid = false;
                cerr << "reciprocity violated\n";
            }

            if (abs(seg->end - seg->endItr->first) > ZERO_CUTOFF)
            {
                // segment end location does not coincide with the first intersection site
                valid = false;
                cerr << "Location of segment end does not coincide with intersection position\n";
            }

            maxDeviation = max(maxDeviation, abs(seg->end - seg->endItr->first));

            if (abs(seg->next()->start - seg->next()->startItr->first) > ZERO_CUTOFF)
            {
                // (next) segment end location does not coincide with the first intersection site
                valid = false;
                cerr << "Location of segment end does not coincide with intersection position\n";
            }
            maxDeviation = max(maxDeviation, abs(seg->next()->start - seg->next()->startItr->first));
        }

        // segment at the end of a trajectory
        else
        {
            // positions should match with wall end location
            if (seg->endItr == seg->trajectory->wallEnd())
            {
                if(abs(seg->end - seg->trajectory->length) > ZERO_CUTOFF)
                {
                    // segment end location does not coincide with the first intersection site
                    valid = false;
                    cerr << "Segment end location does not coincide with end of trajectory\n";
                }
                maxDeviation = max(maxDeviation, abs(seg->end - seg->trajectory->length));
            }

	    // positions should match with wall begin location
            else
            {
                if(abs(seg->end - 0) > ZERO_CUTOFF)
                {
	            // segment end location does not coincide with the first intersection site
                    valid = false;
                    cerr << "Segment end location does not coincide with end of trajectory\n";
                }
                maxDeviation = max(maxDeviation, abs(seg->end - 0));
            }

            if (seg->next()->startItr == seg->next()->trajectory->wallEnd())
            {
                if(abs(seg->next()->start - seg->next()->trajectory->length) > ZERO_CUTOFF)
                {
                    valid = false;
                    cerr << "Segment start location does not coincide with end of trajectory\n";
                    cerr << "start: " << seg->next()->start << ", length: " << seg->next()->trajectory->length << ", diff=" << abs(seg->next()->start - seg->next()->trajectory->length) << "\n";
                }
                maxDeviation = max(maxDeviation, abs(seg->next()->start - seg->next()->trajectory->length));
            }
            else
            {
                if(abs(seg->next()->start - 0) > ZERO_CUTOFF)
                {
                    valid = false;
                    cerr << "Segment start location does not coincide with start of trajectory\n";
                    cerr << "start: " << seg->next()->start << ", length: " << seg->next()->trajectory->length << ", diff=" << abs(seg->next()->start) <<"\n";
                }
                maxDeviation = max(maxDeviation, abs(seg->next()->start - 0));
            }

        }
        seg = seg->next();
    }

    if (!valid)
        cout << "Broken segment connections\n";


    // check begin and end points
    if (plus.dir*plus.velocity > 0)
    {
        if (plus.nextCollision == plus.trajectory->wallEnd())
        {
            // segment end location does not coincide with the first intersection site 
            if (segments.last()->end - plus.trajectory->length > ZERO_CUTOFF)
                valid = false;
        }

        else
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end - plus.nextCollision->first > ZERO_CUTOFF)
                valid = false;
        }

    }

    else
    {
        if (plus.nextCollision == plus.trajectory->wallBegin())
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end < - ZERO_CUTOFF)
                valid = false;
        }
        else
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end - plus.nextCollision->first < - ZERO_CUTOFF)
                valid = false;
        }
    }

    if (!valid)
        cerr << "ERROR: Problem with ordering of plus ends and next collisions.\n";


    // minus end velocity can be zero in the absence of treadmilling, that case skip this check to avoid numerical directionality problems
    if (minus.dir*minus.velocity > ZERO_CUTOFF)
    {
        if (minus.nextCollision == minus.trajectory->wallEnd())
        {
            if (segments.first()->start - minus.trajectory->length > ZERO_CUTOFF)
                valid = false;
        }
        else
        {
            if (segments.first()->start - minus.nextCollision->first > ZERO_CUTOFF)
                valid = false;
        }

    }
    else if (minus.dir*minus.velocity < -ZERO_CUTOFF)
    {
        if (minus.nextCollision == minus.trajectory->wallBegin())
        {
            if (segments.first()->start < - ZERO_CUTOFF)
                valid = false;
        }
        else
        {
            if (segments.first()->start - minus.nextCollision->first < - ZERO_CUTOFF)
                valid = false;
        }
    }

    if (!valid)
        cerr << "ERROR: Problem with ordering of minus ends and next collisions.\n";

    // check validity of intersections by running them to the end of the trajectory
    seg = segments.first();
    while (seg != NULL)
    {
        is = seg->startItr;
        count = 0;
        while ((is != seg->trajectory->wallEnd()) && (count != seg->trajectory->intersections.size()))
        {
            is++;
            count++;
        }
        if (is != seg->trajectory->wallEnd())
        {
            valid = false;
            cerr << "Error: segment start iterator does not belong to correct trajectory\n";
            cerr << "Number of segments in MT: " << segments.size() << "\n";;
        }
        is = seg->endItr;
        count = 0;
        while ((is != seg->trajectory->wallEnd()) && (count != seg->trajectory->intersections.size()))
        {
            is++;
            count++;
        }
        if (is != seg->trajectory->wallEnd())
        {
            valid = false;
            cerr << "Error: segment end iterator does not belong to correct trajectory\n";
            cerr << "Number of segments in MT: " << segments.size() << "\n";;
        }
        seg = seg->next();
    }

    // check the tip's tipmap pointer (round-trip)
    if (system->getEventDescriptor(plus.event.index) != &plus.event)
        valid = false;
    if (system->getEventDescriptor(minus.event.index) != &minus.event)
        valid = false;
    if (system->getEventDescriptor(disappearEvent.index) != &disappearEvent)
        valid = false;

    // round-trip check of the segment-trajectory pairs
    seg = segments.first();
    while (seg != NULL)
    {
        if (*(seg->trajectoryTag) != seg)
            valid = false;

        seg = seg->next();
    }

    // round-trip check of tip-trajectory pairs
    if (*(plus.notificationTag) != &plus)
        valid = false;
    if (*(minus.notificationTag) != &minus)
        valid = false;


    // check orientation of segments / start-end, and tips
    if (segments.first()->dir == minus.dir)
        valid = false;
    if (segments.last()->dir != plus.dir)
        valid = false;

    seg = segments.first();
    while (seg != NULL)
    {
        if (seg->dir == ::forward)
        {
            if (seg->end - seg->start < -ZERO_CUTOFF)
            {
                cerr << "false directionality in positions. diff=" << seg->end - seg->start << "\n";
                valid = false;
            }
            is = seg->startItr;
            count = 0;
            while ((is != seg->endItr) && (count != seg->trajectory->intersections.size()))
            {
                is++;
                count++;
            }
            if (is != seg->endItr)
            {
                cerr << "false directionality in intersection iterators\n";
                valid = false;
            }
        }
        else
        {
            if (seg->end - seg->start > ZERO_CUTOFF)
            {
                cerr << "false directionality in positions. diff=" << seg->end - seg->start << "\n";
                valid = false;
            }
            is = seg->startItr;
            count = 0;
            while ((is != seg->endItr) && (count != seg->trajectory->intersections.size()))
            {
                is--;
                count++;
            }
            if (is != seg->endItr)
            {
                cerr << "false directionality in intersection iterators\n";
                valid = false;
            }
        }

        seg = seg->next();
    }

    // check whether segment end iterators are close to tip->nextCollision

    return valid;
}

Segment::Segment(Microtubule* m, TrajectoryVector& tv)
    : mt(m),
      trajectory(tv.trajectory),
      start(tv.pos),
      end(tv.pos),
      dir(tv.dir),
      trajectoryTag(tv.trajectory->segments.end()),
      nucleationTime(m->system->systemTime + m->system->systemTimeOffset),
      startItr(tv.trajectory->wallEnd()),
      endItr(tv.trajectory->wallEnd())
{
    #ifdef DBG_MTS
    cout << "DBG/MTS: Segment created.\n";
    cout << "begin position: " << start << ", end position: " << end << "\n";
    #endif

    // insert this segment to the associated MT trajectory 
    trajectoryTag = tv.trajectory->insertSegment(this);

    // increase total number of segments for the associated MT
    mt->system->countSegments++;
    return;
}

Segment::~Segment()
{
    #ifdef DBG_MTS
    cout << "DBG/MTS: Segment destroyed.\n";
    #endif

    // decrease total number of segments for the associated MT 
    mt->system->countSegments--;

    // if this is the last segment, remove the tip references...
    trajectory->removeSegment(trajectoryTag);
    return;
}

bool Segment::isLastInMT()
{
    // check whether this is the last segment of the MT
    return (this == mt->segments.last());
}

bool Segment::isFirstInMT()
{
    // check whether this is the first segment of the MT
    return (this == mt->segments.first());
}

bool Segment::crossesIntersection(IntersectionItr& is)
{
    // check whether the segment has a cross intersection with a trajectory
    double temp(0.0);
    temp = (is->first - end) * (is->first - start);

    if (temp < -ZERO_CUTOFF)
        return true;

    else if (temp < ZERO_CUTOFF)
    {
        int sign1 = trajectory->differenceSign((this == mt->segments.last()) ? mt->plus.nextCollision : endItr, end, is, is->first);
        int sign2 = trajectory->differenceSign((this == mt->segments.first()) ? mt->minus.nextCollision : startItr, start, is, is->first);

        if (sign1*sign2 == -1)
            return true;
    }

    return false;
}

MTTip::MTTip(Microtubule* m, TrajectoryVector tv, double v, DeterministicQueue* q, double queueV) :
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
    regionTag = trajectory->base.region->registerOnRegion(this,type(),mt->type);

    // locate the next intersection position of the tip
    locateIntersection();

    // determine the kind of event the tip is going to encounter
    determineEvent();

    return;
}

void MTTip::unlinkFromTrajectory()
{
    // unregister the tip from current region
    trajectory->base.region->unregisterFromRegion(regionTag,type(),mt->type);

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
    if (((dir*velocity) < 0) || ((velocity == 0) && (dir == ::forward)))
        nextCollision--;

    // if direction of the tip is backward then only step down the segment interator
    if (dir == backward)
        temp--;

    // plus tip (must be on the end of top segment)
    if (type() == t_plus)
        mt->segments.last()->endItr = temp;

    // minus tip (must be on the start of bottom segment
    else
        mt->segments.first()->startItr = temp;

    return;
}

void MTTip::advanceIntersection()
{
    // next collission site in front
    if ((dir*velocity > 0) || ((velocity == 0) && (dir == backward)))
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
            temp--;
        else
            temp++;
    }

    // plus tip (must be on the end of top segment)
    if (type() == t_plus)
        mt->segments.last()->endItr = temp;

    // minus tip (must be on the start of bottom segment
    else
        mt->segments.first()->startItr = temp;

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
        else if ((mt->segments.size() != 1) && (this == &(mt->plus)) && (nextCollision == mt->segments.last()->startItr))
        {
            // the segment will be removed
            eventType = ev_end_of_segment;
            eventPos = mt->segments.last()->startItr->first;
        }

        // minus end tip and the next collision is some where at the middle of trajectory
        else if ((mt->segments.size() != 1) && (this == &(mt->minus)) && (nextCollision == mt->segments.first()->endItr))
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
    event.pushOnQueue((eventPos - pos)*dir, eventType);

    // transfer the value of the calculated event position to the next event position 
    nextEventPos = eventPos;

    #ifdef DBG_EVENT
    if (event.type != ev_none)
    {
        cout << "DBG/EVENT: event created. type: " << eventType << ", distance: " << (eventPos - pos)*dir << "\n";
        cout << "tip position: " << position() << ", tip velocity: " << dir*velocity << ", trajectory length: " << trajectory->length << "\n";
        if (nextCollision != trajectory->intersections.end())
            cout << "event position: " << nextCollision->first << "\n";
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
            ref++;
        else
            ref--;
    }

    return;
}

void MTTip::notifyInsert(IntersectionItr& newIs)
{

    // update the MT 
    mt->updateLength();

    IntersectionItr temp = newIs;
    
    // next collission site in front
    if ((dir*velocity > 0) || ((velocity == 0) && (dir==backward)))
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new intersection 
        if ((++temp == nextCollision) && (newIs->first > position()))
        {
            nextCollision = newIs;
            determineEvent();
        }
    }

    // next collission site in back
    else
    {
        // if the scheduled next collision is located after the new intersection point,  then replace it by the new intersection 
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
            refItr = newIs;
    }

    else
    {
        if ((--temp == refItr) && (newIs->first < refPos))
            refItr = newIs;
    }

    return;
}

TipType MTTip::type()
{
    // get tip type of MT
    if (this == &(mt->plus))
        return t_plus;
    else
        return t_minus;
}

double MTTip::position()
{
    // get current tip (top) position on trajectory
    if (this == &(mt->minus))
        return mt->segments.first()->start;
    else
        return mt->segments.last()->end;
}

double MTTip::otherPosition()
{
    // get the (bottom) position of an active segment 
    if (this == &(mt->minus))
        return mt->segments.first()->end;
    else
        return mt->segments.last()->start;
}

Segment& MTTip::segment()
{
    // get the current segement of the tip 
    if (this == &(mt->minus))
        return *(mt->segments.first());
    else
        return *(mt->segments.last());
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
        advanceIntersection();

    // assign an event to the tip
    determineEvent();

    // unregister the tip from old trajectory 
    oldTr->base.region->unregisterFromRegion(regionTag,type(),mt->type);

    // find the new regiontag of the tip
    regionTag = newTr->base.region->registerOnRegion(this,type(),mt->type);

    // copy old notification tag of the tip (required to unregister the tip from old trajectory notification list)
    TrjMTTipTag tempTag = notificationTag;

    // assign new notification tag to the tip
    notificationTag = newTr->registerForNotifications(this);

    // unregister notification of the tip from old trajectory
    oldTr->unregisterForNotifications(tempTag);

    return;
}

