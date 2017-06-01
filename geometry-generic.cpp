/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <list>

#pragma warning(disable:981)
//#pragma warning(disable:177)
#pragma warning(disable:522)
#pragma warning(disable:383)

#include "MersenneTwister.h"
#include "corticalSim.h"


/************************** GENERIC GEOMETRY FUNCTIONS ********************************************/

TrajectoryVector Geometry::createTrajectory(const SurfaceVector& sVec, bool scanExisting)
{
	return sVec.region->insertTrajectory(sVec, scanExisting);
}


TrajectoryVector Geometry::createAndLinkTrajectory(const SurfaceVector& newBase, Trajectory* oldtr, Direction olddir, double cosAngle)
{
	TrajectoryVector tVec = createTrajectory(newBase, true);

	TrajectoryVector thisVec;
	thisVec.trajectory = oldtr;

	if (olddir == ::forward)
	{
		oldtr->nextTr = tVec;
		oldtr->nextTrCosAngle = cosAngle;
		thisVec.dir = backward;
		thisVec.pos = oldtr->length;
	}
	else
	{
		oldtr->prevTr = tVec;
		oldtr->prevTrCosAngle = cosAngle;
		thisVec.dir = ::forward;
		thisVec.pos = 0;
	}

	if (tVec.dir == ::forward)
	{
		tVec.trajectory->prevTr = thisVec;
		tVec.trajectory->prevTrCosAngle = cosAngle;
	}
	else
	{
		tVec.trajectory->nextTr = thisVec;
		tVec.trajectory->nextTrCosAngle = cosAngle;
	}

	return tVec;
}


SurfaceVector Geometry::randomSurfaceVector()
{
	int regionSelector = 0;

	if (regions.size() > 1)
	{
		double rArea = system->randomGen.randDblExc(area);
	
		if (rArea < 0.5*area)
		{
			regionSelector = 0;
			while (rArea > regions[regionSelector]->area)
			{
				rArea -= regions[regionSelector++]->area;
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
			regionSelector = regions.size()-1;
			rArea = area - rArea;
			while (rArea > regions[regionSelector]->area)
			{
				rArea -= regions[regionSelector--]->area;
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
	
	return regions[regionSelector]->randomSurfaceVector();
}


bool Geometry::integrityCheck()
{
	int ridx = 0;
	Trajectory* tptr;
	bool valid = true;
	
	int plusGrowCount = 0;
	int plusShrinkCount = 0;
	int minusCount = 0;

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

	if (plusGrowCount != system->growingTipsNormal + system->growingTipsSpecial)
	{
		cerr << "Integrity check failed: number of growing tips does not equal cached values\n";
		valid = false;
	}

	return valid;	
}

double Geometry::opticalLength()
{
	int ridx = 0;
	double sum = 0.0;
	
	for (ridx = 0; ridx < regions.size(); ridx++)
	{
		sum += regions[ridx]->opticalLength();
	}
	return sum;	
}

int Geometry::trajectoryCount()
{
	int ridx = 0;
	int sum = 0;
	
	for (ridx = 0; ridx < regions.size(); ridx++)
	{
		sum += regions[ridx]->trajectories.size();
	}
	return sum;	
}

/************************* GENERIC REGION FUNCTIONS **************************************************/

TrajectoryVector Region::insertTrajectory(const SurfaceVector& sVec, bool scanExisting)
{
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Region::insertTrajectory() called\n";
	#endif
	TrajectoryVector tVec;
	double tLength;

	SurfaceVector newBase = sVec;
	getTrajectoryCoordinates(newBase, tLength, tVec);

	if (scanExisting)
	{
		Trajectory* tr = trajectories.first();
		while (tr != NULL)
		{
			// Region::getTrajectoryCoordinates() always generates the same surface vector coordinates for inputs on the same physical trajectory.
			// The orientation of the trajectory itself is therefore also fixed.
			if ((abs(tr->base.x - newBase.x) < ZERO_CUTOFF) && (abs(tr->base.y - newBase.y) < ZERO_CUTOFF) && (abs(tr->base.angle - newBase.angle) < ZERO_CUTOFF))
			{	
				//a suitable trajectory already exists. Now check if it's free to connect
				if (((tVec.dir == ::forward) && (tr->prevTr.trajectory == NULL)) ||
					((tVec.dir == ::backward) && (tr->nextTr.trajectory == NULL)))
				{
#ifdef DBG_ASSERT
					cerr << "Region::insertTrajectory: Existing trajectory found.\n";
#endif	
					tVec.trajectory = tr;
					return tVec;
				}
#ifdef DBG_ASSERT
					cerr << "Region::insertTrajectory: Existing trajectory found, but already connected. [Likely a bug]\n";
#endif	
			}

			tr = tr->next();
		}
	}

	tVec.trajectory = trajectories.create(newBase,tLength);
	makeIntersectionList(tVec.trajectory);

	geometry->system->countIntersections += 2*(tVec.trajectory->intersections.size());
	return tVec;
}



void Region::removeTrajectory(Trajectory* tr)
{
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Region::removeTrajectory() called\n";
	#endif
	
	#ifdef DBG_ASSERT
	if (!tr->segments.empty())
	{
		cout << "Big problem";
	}
	if (!tr->notificationList.empty())
	{
		cout << "ojee...\n";
	}
	#endif

	geometry->system->countIntersections -= 2* (tr->intersections.size());

	Trajectory* otherTr;
	IntersectionItr is;
	
	is = ++tr->intersections.begin();	//start one past the beginning since the first element is a stub...
	while (is != tr->intersections.end())
	{
		otherTr = is->second.otherTrajectory;
		otherTr->invalidateIntersection(is->second.mirror);
		otherTr->intersections.erase(is->second.mirror);
		is++;
	}
	
	tr->intersections.clear();


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

	trajectories.remove(tr);
}

double Region::opticalLength()
{
	Trajectory* tr;
	double sum = 0.0;
	
	tr = trajectories.first();
	while (tr != NULL)
	{
		sum += tr->coveredLength();
		tr = tr->next();
	}
	return sum;	
}



RegionMTTipTag Region::registerOnRegion(MTTip* pTip, TipType tiptype, MTType mttype)
{
	updateRegionLength();

	if (tiptype == t_minus)
		return minusTipList.insert(minusTipList.end(),pTip);
	else
	{
		if (mttype == mt_growing)
		{
			if (parameterType == r_param_normal)
				geometry->system->growingTipsNormal++;
			else // if parameterType == r_param_special
				geometry->system->growingTipsSpecial++;
			return growingPlusTipList.insert(growingPlusTipList.end(),pTip);
		}
		else
			return shrinkingPlusTipList.insert(shrinkingPlusTipList.end(),pTip);
	}
}

void Region::unregisterFromRegion(RegionMTTipTag tag, TipType tiptype, MTType mttype)
{
	updateRegionLength();

	if (tiptype == t_minus)
		minusTipList.erase(tag);
	else
	{
		if (mttype == mt_growing)
		{
			if (parameterType == r_param_normal)
				geometry->system->growingTipsNormal--;
			else // if parameterType == r_param_special
				geometry->system->growingTipsSpecial--;
			growingPlusTipList.erase(tag);
		}
		else
			shrinkingPlusTipList.erase(tag);
	}
	return;
}


/****************************** TRAJECTORY FUNCTIONS **************************/

Trajectory::Trajectory(SurfaceVector baseVec, double l) : 
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
	base.region->geometry->system->countTrajectories++;
	intersections.insert(pair<double, Intersection>(-1.0, Intersection()));

	return;	
}

Trajectory::~Trajectory()
{
	base.region->geometry->system->countTrajectories--;
	return;
}

bool Trajectory::integrityCheck()
{
	bool valid = true;
	
	list<Segment*>::iterator segItr;
	IntersectionItr is;

	int isCount;	
	
	
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
	is = wallBegin();
	for (isCount = 1; isCount < intersections.size(); isCount++)
	{
		is++;
		if (is->second.mirror->second.mirror != is)
		{
			cout << "Integrity: intersection mirroring broken.\n";
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
	double sumLength = 0.0;
	list<Segment*>::iterator seg;
	
	seg = segments.begin();
	while (seg != segments.end())
	{
		sumLength += (**seg).length();
		seg++;
	}
	return sumLength;
}


double Trajectory::coveredLength()
{
	list<double> startPoints;
	list<double> endPoints;
	list<double>::iterator stepup;
	list<double>::iterator stepdown;
	
	list<Segment*>::iterator seg;
	
	seg = segments.begin();
	while (seg != segments.end())
	{
		startPoints.push_back(min((**seg).start, (**seg).end));
		endPoints.push_back(max((**seg).start, (**seg).end));
		seg++;
	}
	startPoints.push_back(VERY_LARGE);
	startPoints.sort();
	endPoints.sort();

	stepup = startPoints.begin();
	stepdown = endPoints.begin();
	double previousPos = 0;
	double currentPos;
	int occupancy = 0;
	double cLength=0;

	while (stepdown != endPoints.end())
	{
		currentPos = min(*stepup, *stepdown);
		if (occupancy >= 1)
			cLength += currentPos - previousPos;
		previousPos = currentPos;
		if ((*stepup) < (*stepdown))
		{
			stepup++;
			occupancy++;
		}
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
	list<MTTip*>::iterator tip;

	// notify all tips on the current trajectory of invalidation	
	tip = notificationList.begin();
	while (tip != notificationList.end())
	{
		(**tip).notifyRemove(oldIs);
		tip++;
	}

	return;
}


void Trajectory::newIntersection(IntersectionItr& newIs)
/* called from makeIntersectionList()
 * 
 */
{
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Trajectory::newIntersection() called.\n";
	#endif

	list<MTTip*>::iterator tip;

	IntersectionItr compare;
	
	int occupancy = 0;
	
	list<Segment*>::iterator seg;
	seg = segments.begin();
	while (seg != segments.end())
	{
		(**seg).mt->updateLength();

		if ((**seg).crossesIntersection(newIs))
		{	occupancy++;	}

		seg++;
	}
	
	newIs->second.occupancy = occupancy;

	// notify all tips on the current trajectory of insertion
	tip = notificationList.begin();
	while (tip != notificationList.end())
	{
		(**tip).notifyInsert(newIs);
		tip++;
	}


	return;
}

int Trajectory::differenceSign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2)
/* returns the sign of 'pos1 - pos2'. If the two positions are equal, the iterators are used instead to determine order.
 */
{
	if (pos1 - pos2 > ZERO_CUTOFF)
		return 1;
	else if (pos1 - pos2 < -ZERO_CUTOFF)
		return -1;
	else
	{
		if (itr1 == itr2)
			return 0;
		if (itr1 == wallEnd())
			return 1;

		do
		{
			itr2++;
		}
		while ((itr2 != wallEnd()) && (itr1 != itr2));
		
		if (itr1 == itr2)
			return 1;
		else
			return -1;
	}	
}

TrjSegmentTag Trajectory::insertSegment(Segment* s)
{
	#ifdef DBG_ASSERT
	// check whether the size of the inserted segment equals zero
	// for insertion of non-zero segments, care must be taken to update the occupancy numbers
	if (abs(s->end - s->start) > ZERO_CUTOFF)
	{
		cerr << "DBG/ASSERT: ERROR: Not permitted to insert segment with non-zero length.\n";
		exit(-1);
	}
	#endif
	return segments.insert(segments.end(),s);
}

void Trajectory::removeSegment(TrjSegmentTag s)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Trajectory::removeSegment() called.\n";
#endif

	segments.erase(s);

	conditionalRemove();
	return;	
}


TrjMTTipTag Trajectory::registerForNotifications(MTTip* pTip)
{
	return notificationList.insert(notificationList.end(),pTip);
}

void Trajectory::unregisterForNotifications(TrjMTTipTag tag)
{
	notificationList.erase(tag);
	conditionalRemove();
	return;
}




TrajectoryVector Trajectory::nextTrajectory(Direction dir)
{
	if (dir == ::forward)
	{
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
		if (prevTr.trajectory != NULL)
		{
			#ifdef DBG_GEOMETRY
			cout << "DBG/GEOMETRY: Following existing trajectory link forward. [from " << RegionTypeText[base.region->type] << \
				" to " << RegionTypeText[prevTr.trajectory->base.region->type] << "]\n";
			
			#endif
			return prevTr;
		}
	}
	
	
	return base.region->geometry->extendTrajectory(this,dir);
}
