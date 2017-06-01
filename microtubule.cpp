/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <cmath>
#include <cstdlib>


#pragma warning(disable:981)
//#pragma warning(disable:177)
#pragma warning(disable:522)
#pragma warning(disable:383)

#include "MersenneTwister.h"
#include "corticalSim.h"



/******************************** Microtubule functions *************************/

Microtubule::Microtubule(System* s, TrajectoryVector tv, bool initialize) 
		: system(s),
		  plus(this, tv, s->p.vPlus, &(s->vPlusQueue), s->p.vPlus),
		  minus(this, tv.flipped(), -s->p.vTM, &(s->timeQueue), -s->p.vTM),		//NOTE: REVERSE DIRECTION
		  disappearEvent(this, &(s->timeQueue), - s->p.vMin + s->p.vTM),
		  type(mt_growing),
		  nucleationTime(s->systemTime + s->systemTimeOffset),
		  previousUpdateTag(s->currentTimeTag)
/* NOTE: some compilers may complain about the presence of 'this' in the initializer
* list. This is no problem as long as the 'this' pointer is not used to access any of 
* the members during the construction of the object
*/
{
	#ifdef DBG_MTS
	cout << "DBG/MTS: Microtubule created.\n";
	#endif

	if (initialize)
	{
		segments.create(this,tv);
		plus.initialize();
		minus.initialize();
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
	
	system->mtLifetimeHistogram.insert(system->systemTime + system->systemTimeOffset - nucleationTime);

	plus.unlinkFromTrajectory();
	minus.unlinkFromTrajectory();

	return;	
}


void Microtubule::setDisappearEvent()
{
	if ((segments.size() == 1) && (type == mt_shrinking)) // && (plus.velocity - minus.velocity < 0)
	{
		disappearEvent.pushOnQueue(segments.first()->length(), ev_disappear);
	}
	else
		disappearEvent.clear();
	return;	
}


void Microtubule::handleEvent(const EventDescriptor* ei) 
{

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
	case ev_wall:
		wall();
		break;
	case ev_collision:
		collision();
		break;
	case ev_backtrack:
		if (ei == &(plus.event))
			backtrack(&plus);
		else
			backtrack(&minus);
		break;
	case ev_end_of_segment:
		if (ei == &(plus.event))
			endOfSegment(&plus);
		else
			endOfSegment(&minus);
		break;
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
	
	updateLength();
	
	plus.trajectory->base.region->unregisterFromRegion(plus.regionTag,t_plus,mt_growing);
	type = mt_shrinking;
	plus.regionTag = plus.trajectory->base.region->registerOnRegion(&plus,t_plus,mt_shrinking);
	plus.velocity = system->p.vMin;
	plus.event.reinitialize(&(system->timeQueue),system->p.vMin);
	system->shrinking_mts.import(system->growing_mts, this);

	if (system->p.vMin < 0)
	{
		// velocity has changed sign... need to select a new event
		plus.advanceIntersection();
	}
	
	// OPTIONAL: if this is triggered by a collision, and vMin > 0, increase occupancy
	
	plus.determineEvent();
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
	

	plus.trajectory->base.region->unregisterFromRegion(plus.regionTag,t_plus,mt_shrinking);
	type = mt_growing;
	plus.regionTag = plus.trajectory->base.region->registerOnRegion(&plus,t_plus,mt_growing);
	plus.velocity = system->p.vPlus;
	plus.event.reinitialize(&(system->vPlusQueue), system->p.vPlus);
	system->growing_mts.import(system->shrinking_mts, this);

	if (system->p.vMin < 0)
	{
		// velocity has changed sign... need to select a new event
		plus.advanceIntersection();
	}
	plus.determineEvent();
	setDisappearEvent();

	return;
}

void Microtubule::sever(Segment* cutSeg, double cutPos)
{
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
	double splitPos = is->first;
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
	Microtubule* newMT;
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
// only called for plus end
{
	TrajectoryVector tv;
	
	// set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
	segments.last()->end = plus.nextEventPos;
	
	
	tv = plus.trajectory->nextTrajectory(plus.dir);

	// basic 'forbidden zone' catastrophes
	if ((system->p.forbiddenZones) && (tv.trajectory->base.region->parameterType == r_param_modified))
	{
		tv.trajectory->conditionalRemove();
		catastrophe();
		return;
	}

	// edge catastrophe catastrophes
	if (system->p.edgeCatastropheEnabled)
	{
		double cosAngle = (plus.dir == ::forward) ? plus.trajectory->nextTrCosAngle : plus.trajectory->prevTrCosAngle;

		//NOTE: trajectories that traverse from 'normal' to 'special' regions without affecting the orientation (cosAngle==1) will NOT induce catastrophes.
		if (1.0 - cosAngle > ZERO_CUTOFF)
		{
			// if we transition from a 'normal' to 'special' region or vice versa, use pCatSpecialEdge
			double pCat;
			if (((plus.trajectory->base.region->parameterType == r_param_normal) && (tv.trajectory->base.region->parameterType == r_param_modified))
				|| ((plus.trajectory->base.region->parameterType == r_param_modified) && (tv.trajectory->base.region->parameterType == r_param_normal)))
				pCat = system->p.pCatSpecialEdge;
			else
				pCat = system->p.pCatRegularEdge;

			// is smooth?
			if (system->p.edgeCatastropheSmooth)
				pCat *= (1.0 - pow(cosAngle,2));

			if (system->randomGen.rand() <= pCat)
			{
				tv.trajectory->conditionalRemove();
				catastrophe();
				return;
			}


		}

		
	}

	system->boundaryCrossingCount++;


	#ifdef DBG_ASSERT
	if 	(segments.last()->endItr != plus.nextCollision)
		cerr << "DBG/ASSERT: nextCollision does not match with endItr upon entering wall event.\n";
	#endif
	segments.last()->endItr = plus.nextCollision;
	segments.create(this, tv);
	segments.last()->startItr = tv.dir == ::forward ? tv.trajectory->wallBegin() : tv.trajectory->wallEnd();
	plus.switchTrajectory(tv.trajectory, tv.dir, tv.dir == ::forward ? tv.trajectory->wallBegin() : tv.trajectory->wallEnd());
	
	if (segments.size() == 2)
	{
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

	if (plus.nextCollision->second.mirror->second.occupancy == 0)
		type = ct_crossover;
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
	
		type = system->collisionResult(diffAngle, plus.nextCollision->second.occupancy, plus.nextCollision->second.mirror->second.occupancy);
	}

	switch(type)
	{
	case ct_zipper: 
		system->totalZipperCount++;
		zipper(dir);
		break;
	case ct_crossover:
		system->totalCrossoverCount++;
		crossover();
		break;
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

	TrajectoryVector tv(plus.nextCollision->second.mirror->first, dir, plus.nextCollision->second.otherTrajectory);

	segments.last()->endItr = plus.nextCollision;
	segments.create(this,tv);
	segments.last()->startItr = plus.nextCollision->second.mirror;
	plus.switchTrajectory(plus.nextCollision->second.otherTrajectory, dir, plus.nextCollision->second.mirror);

	if (segments.size() == 2)
	{
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
		exit(-1);
	}
	#endif
	
	Direction newDir;
	Segment* newSeg;
	Segment* killSeg;

	if (tip == &plus)
	{
		killSeg = segments.last();
		newSeg = segments.last()->previous();
		newDir = newSeg->dir;
	}
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

	segments.remove(killSeg);


	if ((tip->nextCollision == tip->trajectory->wallBegin()) ||(tip->nextCollision == tip->trajectory->wallEnd()))
	{
		tip->switchTrajectory(newSeg->trajectory, newDir, newDir == ::forward ? newSeg->trajectory->wallEnd() : newSeg->trajectory->wallBegin());
		system->boundaryCrossingCount--;
	}
	else
		tip->switchTrajectory(tip->nextCollision->second.otherTrajectory, newDir, tip->nextCollision->second.mirror);

	if (segments.size()==1)
	{
		if (tip == &(plus))
			minus.determineEvent();
		else
			plus.determineEvent();
		setDisappearEvent();
	}		

	return;
}

// this may become a tip member function
void Microtubule::backtrack(MTTip* tip)
{
	if (tip == &(plus))
		segments.last()->end = plus.nextEventPos;
	else
		segments.first()->start = minus.nextEventPos;

	tip->nextCollision->second.occupancy--;
	if ((tip->nextCollision->second.occupancy == 0) && (tip->nextCollision->second.mirror->second.occupancy > 0))
	{
		system->removeOccupiedIntersection(tip->nextCollision->second);
	}

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
	tip->advanceIntersection();
	tip->determineEvent();
	return;
}

// this may become a tip member function
void Microtubule::crossover()
{
	#ifdef DBG_ASSERT
	if (plus.nextCollision->second.occupancy < 0)
	{
		cerr << "DBG/ASSERT: ERROR: negative intersection occupancy before crossover.\n";
		exit(-1);
	}
	#endif

	if ((plus.nextCollision->second.occupancy == 0) && plus.nextCollision->second.mirror->second.occupancy > 0)
	{
		system->addOccupiedIntersection(plus.nextCollision);
	}
	plus.nextCollision->second.occupancy++;
	plus.advanceIntersection();
	plus.determineEvent();
	return;
}

void Microtubule::harakiri()
{
	if (type == mt_growing)
		system->growing_mts.remove(this);
	else
		system->shrinking_mts.remove(this);

	return;
}

bool Microtubule::integrityCheck()
{
	bool valid = true;
	double maxDeviation = 0.0;
	
	Segment* seg;
	IntersectionItr is;
	int count;
	
	// check connections between end-of segment iterators, and match with positions (should be exact)

	if (segments.size() == 0)
	{
		valid = false;
		cerr << "Microtubule contains no segments.\n";
	}
	else if (segments.size() != 1)
	{
		if (segments.size() == 0)
		{
			cerr  << "microtubule without segments!\n";
			exit(-1);
		}
		seg = segments.first();
		if (seg == NULL)
		{
			cerr << "null pointer received\n";
			exit(-1);
		}
		
		if ((seg->endItr != seg->trajectory->wallEnd()) && (seg->endItr != seg->trajectory->wallBegin()))
		{
			if (seg->endItr->second.mirror != seg->next()->startItr)
			{
				valid = false;
				cerr << "reciprocity violated\n";
			}
			if (abs(seg->end - seg->endItr->first) > ZERO_CUTOFF)
			{
				valid = false;
				cerr << "Location of segment end does not coincide with intersection position\n";
			}
			maxDeviation = max(maxDeviation, abs(seg->end - seg->endItr->first));
			if (abs(seg->next()->start - seg->next()->startItr->first) > ZERO_CUTOFF)
			{
				valid = false;
				cerr << "Location of segment end does not coincide with intersection position\n";
			}
			maxDeviation = max(maxDeviation, abs(seg->next()->start - seg->next()->startItr->first));
		}
		else
		{
			// positions should match with wall
			if (seg->endItr == seg->trajectory->wallEnd())
			{
				if(abs(seg->end - seg->trajectory->length) > ZERO_CUTOFF)
				{
					valid = false;
					cerr << "Segment end location does not coincide with end of trajectory\n";
				}
				maxDeviation = max(maxDeviation, abs(seg->end - seg->trajectory->length));
			}
			else
			{
				if(abs(seg->end - 0) > ZERO_CUTOFF)
				{
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
			if (segments.last()->end - plus.trajectory->length > ZERO_CUTOFF)
				valid = false;
		}
		else
		{
			if (segments.last()->end - plus.nextCollision->first > ZERO_CUTOFF)
				valid = false;	
		}
		
	}
	else
	{
		if (plus.nextCollision == plus.trajectory->wallBegin())
		{
			if (segments.last()->end < - ZERO_CUTOFF)
				valid = false;
		}
		else
		{
			if (segments.last()->end - plus.nextCollision->first < - ZERO_CUTOFF)
				valid = false;	
		}
	}
	if (!valid)
		cerr << "ERROR: Problem with ordering of plus ends and next collisions.\n";

	// check begin and end points
	// Note that the minus end velocity can be zero in the absence of treadmilling.
	// In that case: skip this check to avoid numerical directionality problems.
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

//	cout << "Check performed: maximum deviation is " << maxDeviation << "\n";
	
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


	// 
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
	
	// suggestion for future expansion: check whether segment end iterators are close to tip->nextCollision
	
	return valid;
}


/******************************** Segment functions ***************************/


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
	trajectoryTag = tv.trajectory->insertSegment(this);
	mt->system->countSegments++;
	return;
}


Segment::~Segment()
{
	#ifdef DBG_MTS
	cout << "DBG/MTS: Segment destroyed.\n";
	#endif

	mt->system->countSegments--;
	// store segment lifetime
	if (trajectory->base.region->accountingType == r_accounting_normal)		
		mt->system->segAngleLifetimeHistogram.insert(trajectory->base.angle, mt->system->systemTime + mt->system->systemTimeOffset - nucleationTime);

	trajectory->removeSegment(trajectoryTag);
	return;
}

bool Segment::isLastInMT()
{ 
	return (this == mt->segments.last()); 
}

bool Segment::isFirstInMT()
{ 
	return (this == mt->segments.first()); 
}

bool Segment::crossesIntersection(IntersectionItr& is)
{
	double temp;

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
/*************************** MT tip functions ******************************/



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
	notificationTag = trajectory->registerForNotifications(this);
	regionTag = trajectory->base.region->registerOnRegion(this,type(),mt->type);
	locateIntersection();
	determineEvent();
	return;
}

void MTTip::unlinkFromTrajectory()
{
	trajectory->base.region->unregisterFromRegion(regionTag,type(),mt->type);
	//this one should come last, because it may lead to the removal of the trajectory.
	trajectory->unregisterForNotifications(notificationTag);
	return;
}


void MTTip::locateIntersection()
{
	// go through the intersection list looking for the first event that can take place
	#ifdef DBG_MTS
	cout << "DBG/MTS: MTTip::locateIntersection()\n";
	#endif
	
	// find the first intersection after 'position()'
	nextCollision = trajectory->intersections.upper_bound(position());
	IntersectionItr temp = nextCollision;
	// if necessary, step to the intersection on the other side of the tip
	if (((dir*velocity) < 0) || ((velocity == 0) && (dir == ::forward)))
		nextCollision--;

	if (dir == backward)
		temp--;

	if (type() == t_plus)
		mt->segments.last()->endItr = temp;
	else
		mt->segments.first()->startItr = temp;

	return;
}


void MTTip::advanceIntersection()
{
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
	if (velocity <= 0)
	{
		if (dir == backward)
			temp--;
		else
			temp++;
	}

	if (type() == t_plus)
		mt->segments.last()->endItr = temp;
	else
		mt->segments.first()->startItr = temp;
	
	return;
}


void MTTip::determineEvent()
{
	double eventPos;
	DeterministicEventType eventType;
	double pos = position();
	
	if (velocity > 0)
	{
		if (nextCollision == trajectory->wallBegin())
		{
			eventType = ev_wall;
			eventPos = 0;
		}
		else if (nextCollision == trajectory->wallEnd())
		{
			eventType = ev_wall;
			eventPos = trajectory->length;
		}
		else
		{
			eventType = ev_collision;
			eventPos = nextCollision->first;
		}
	}
	else // (velocity<=0), i.e. a tip that's 'eating' tubulin
	{
		if (nextCollision == trajectory->wallBegin())
		{
			eventType = ev_end_of_segment;
			eventPos = 0;
		}
		else if (nextCollision == trajectory->wallEnd())
		{
			eventType = ev_end_of_segment;
			eventPos = trajectory->length;
		}
		else if ((mt->segments.size() != 1) && (this == &(mt->plus)) && (nextCollision == mt->segments.last()->startItr))
		{
			eventType = ev_end_of_segment;
			eventPos = mt->segments.last()->startItr->first;
		}
		else if ((mt->segments.size() != 1) && (this == &(mt->minus)) && (nextCollision == mt->segments.first()->endItr))
		{
			eventType = ev_end_of_segment;
			eventPos = mt->segments.first()->endItr->first;
		}
		else
		{
			eventType = ev_backtrack;
			eventPos = nextCollision->first;
		}
	}


	event.pushOnQueue((eventPos - pos)*dir, eventType);
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
	mt->updateLength();
	
	if (oldIs == nextCollision)
	{
		advanceIntersection();
		determineEvent();
	}
	
	IntersectionItr& ref = (type() == t_plus) ? mt->segments.last()->endItr : mt->segments.first()->startItr;
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
	
	mt->updateLength();

	IntersectionItr temp = newIs;
	if ((dir*velocity > 0) || ((velocity == 0) && (dir==backward)))
	{
		if ((++temp == nextCollision) && (newIs->first > position()))
		{
			nextCollision = newIs;
			determineEvent();	
		}
	}
	else
	{
		if ((--temp == nextCollision) && (newIs->first < position()))
		{
			nextCollision = newIs;
			determineEvent();	
		}
	}

	temp = newIs;
	IntersectionItr& refItr = (type() == t_plus) ? mt->segments.last()->endItr : mt->segments.first()->startItr;
	double& refPos = (type() == t_plus) ? mt->segments.last()->end : mt->segments.first()->start;
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
	if (this == &(mt->plus)) 
		return t_plus; 
	else 
		return t_minus; 
}

double MTTip::position()
{
	if (this == &(mt->minus))
		return mt->segments.first()->start;
	else 
		return mt->segments.last()->end;
}

double MTTip::otherPosition()
{
	if (this == &(mt->minus))
		return mt->segments.first()->end;
	else 
		return mt->segments.last()->start;
}

Segment& MTTip::segment()
{
	if (this == &(mt->minus))
		return *(mt->segments.first());
	else
		return *(mt->segments.last());
}

void MTTip::switchTrajectory(Trajectory* newTr, Direction d, IntersectionItr intersect, bool advance)
{
	Trajectory* oldTr = trajectory;

	trajectory = newTr;
	dir = d;
	nextCollision = intersect;

	if (advance)
		advanceIntersection();
	determineEvent();

	oldTr->base.region->unregisterFromRegion(regionTag,type(),mt->type);
	regionTag = newTr->base.region->registerOnRegion(this,type(),mt->type);

	TrjMTTipTag tempTag = notificationTag;
	notificationTag = newTr->registerForNotifications(this);
	// call the following functions *after* the previous statements, because it can cause
	// the trajectory - and therefore the intersection array - to be deleted. This invalidates the 
	// initial 'nextCollision' iterator and the end of segment iterator
	oldTr->unregisterForNotifications(tempTag);
	return;
}

