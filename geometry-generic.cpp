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


Aster Geometry::createAster(int nr, SurfaceVector sv) 
{	
	return sv.region->createAster(nr,sv);
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


void Geometry::removeMetaTrajectory(MetaTrajectory* mtr)
{
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Geometry::removeMetaTrajectory() called\n";
	#endif
	metaTrajectories.remove(mtr);
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
		otherTr = is->ISREF.otherTrajectory;
		otherTr->invalidateIntersection(is->ISREF.mirror);
		otherTr->intersections.erase(is->ISREF.mirror);
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
		if (is->ISREF.mirror->ISREF.mirror != is)
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
				is->ISREF.occupancy--;
		}
		else
		{
			while (--is != (**segItr).endItr)
				is->ISREF.occupancy--;
		}
		segItr++;
	}

	// now, check whether all occupancy counts are zero
	is = intersections.begin();
	while (++is != intersections.end())
	{
		if (is->ISREF.occupancy != 0)
		{
			cerr << "Detected occupancy " << is->ISREF.occupancy << " (should be zero).\n";
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
				is->ISREF.occupancy++;
		}
		else
		{
			while (--is != (**segItr).endItr)
				is->ISREF.occupancy++;
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

// Working version
/*
void Trajectory::coveredLengthSplit(double * length, double* Olength, double* cuts, int n,int rev)
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
	double cOLength=0;
	int curL,shift;

	if (rev > 0)
	{
		curL=0;
		shift=1;
	}
	else
	{
		curL=n-1;
		shift = 0;
	}
	while (cuts[curL+shift] < 0)
		curL += rev;
 
	// TODO: TO COUNT NUMBER OF MTs: make pairs: startPoint/endPoint + weight = 1/mt.length. (in stead of unifor weight 1 for length. 

	length[0] = 0.;
	while (stepdown != endPoints.end())
	{
		currentPos = min(*stepup, *stepdown);
		while (curL >= 0 && curL < n+shift && cuts[curL+shift] < previousPos) {
			curL += rev;
			if (curL < 0) {
				cout << "Should not happen! " << rev << " " << cuts[1] << " " << curL << "\n";
				curL=0;
				break;
			}
		}
		while (currentPos > cuts[curL+shift] ) {
			if (occupancy >= 1) {
				cLength += (cuts[curL+shift] - previousPos) * occupancy;
				cOLength += cuts[curL+shift]  - previousPos;
				previousPos = cuts[curL+shift];
				Olength[curL] = cOLength;
				length[curL] = cLength;
				cLength = 0.;
				cOLength = 0.;
				//cLength = (min(currentPos,cuts[curL+shift+rev]) - previousPos) * occupancy;
				//cOLength = min(currentPos,cuts[curL+shift+rev]) - previousPos;
			}		
			else {
				length[curL] = cLength;
				cLength = 0.;
				Olength[curL] = cOLength;
				cOLength = 0.;
			}
			curL += rev;
			if (curL < 0)
			{
				cout << "Should not happen either! " << rev << " " << cuts[0] << " " << curL << "\n";
				curL = 0;
			}
			if (curL > n) {
				cerr << "function Trajectory::coveredLengthSplit(...): out of bounds! n=" << n << " < " << curL << "\n";
				for (int k=0 ; k<n+1; k++) {
					cerr << cuts[k] << "\n";
				}
				exit(35);
			}
		}
		//else
		if (occupancy >= 1)
		{
			cLength += (currentPos - previousPos) * occupancy;
			cOLength += currentPos - previousPos;
		}
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
#ifdef DBG_SPATIAL
	if (currentPos > cuts[curL+shift])
	{
		cout << "lost something " << previousPos << " " << currentPos << " " << cuts[curL+shift] << " " << curL << " " << shift << " " << rev <<  " " << n << "\n";
		for (int k=0 ; k<n+1; k++) {
			cerr << cuts[k] << "\n";
		}
	}
#endif
	length[curL] = cLength;
	Olength[curL] = cOLength;
	return;
}
*/

//* Broken version; aiming for new feature.
void Trajectory::coveredLengthSplit(double * length, double* Olength, double* count, double* cuts, int n,int rev)
{
	list<pair <double,double>> startPoints;
	list<pair <double,double>> endPoints;
	list<pair <double,double>>::iterator stepup;
	list<pair <double,double>>::iterator stepdown;
	double weight;
	
	list<Segment*>::iterator seg;
	
	seg = segments.begin();
	while (seg != segments.end())
	{
		//weight = 1./((**seg).mt->length()+ZERO_CUTOFF); // gives a small roundoff error, but prevents zero division errors for just nucleated MTs. 
		weight = 1./((**seg).mt->length()); // gives a small roundoff error, but prevents zero division errors for just nucleated MTs. 
		startPoints.push_back(make_pair(min((**seg).start, (**seg).end),weight));
		endPoints.push_back(make_pair(max((**seg).start, (**seg).end),weight));
		seg++;
	}
	startPoints.push_back(make_pair(VERY_LARGE,0.));
	startPoints.sort();
	endPoints.sort();

	stepup = startPoints.begin();
	stepdown = endPoints.begin();
	double previousPos = 0;
	double currentPos;
	double currentWeight = 0.;
	int occupancy = 0;
	double cLength=0;
	double cCount=0;
	double cOLength=0;
	int curL,shift;

	if (rev > 0)
	{
		curL=0;
		shift=1;
	}
	else
	{
		curL=n-1;
		shift = 0;
	}
	while (cuts[curL+shift] < 0)
		curL += rev;
 
	// TODO: TO COUNT NUMBER OF MTs: make pairs: startPoint/endPoint + weight = 1/mt.length. (in stead of unifor weight 1 for length. 

	length[0] = 0.;
	while (stepdown != endPoints.end())
	{
		currentPos = min((*stepup).first, (*stepdown).first);
		while (curL >= 0 && curL < n+shift && cuts[curL+shift] < previousPos) {
			curL += rev;
			if (curL < 0) {
				cout << "Should not happen! " << rev << " " << cuts[1] << " " << curL << "\n";
				curL=0;
				break;
			}
		}
		while (currentPos > cuts[curL+shift] ) {
			if (occupancy >= 1) {
				cLength += (cuts[curL+shift] - previousPos) * occupancy;
				cOLength += cuts[curL+shift]  - previousPos;
				cCount += (cuts[curL+shift] - previousPos) * currentWeight;
				previousPos = cuts[curL+shift];
				Olength[curL] += cOLength;
				length[curL] += cLength;
				count[curL] += cCount;
				cLength = 0.;
				cOLength = 0.;
				cCount = 0.;
				//cLength = (min(currentPos,cuts[curL+shift+rev]) - previousPos) * occupancy;
				//cOLength = min(currentPos,cuts[curL+shift+rev]) - previousPos;
			}		
			else {
				length[curL] += cLength;
				cLength = 0.;
				Olength[curL] += cOLength;
				cOLength = 0.;
				count[curL] += cCount;
				cCount = 0.;
			}
			curL += rev;
			if (curL < 0)
			{
				cout << "Should not happen either! " << rev << " " << cuts[0] << " " << curL << " " << currentPos << "\n";
				for (int k=0 ; k<n+1; k++) {
					cerr << cuts[k] << "\n";
				}
				curL = 0;
			}
			if (curL > n) {
				cerr << "function Trajectory::coveredLengthSplit(...): out of bounds! n=" << n << " < " << curL << "\t"<< rev << "\t"<< cuts[curL+shift] << "\t" << currentPos << "\n";
				for (int k=0 ; k<n+1; k++) {
					cerr << cuts[k] << "\n";
				}
				exit(35);
			}
		}
		//else
		if (occupancy >= 1)
		{
			cLength += (currentPos - previousPos) * occupancy;
			cCount += (currentPos - previousPos) * currentWeight;
			cOLength += currentPos - previousPos;
		}
		previousPos = currentPos;
		if ((*stepup).first < (*stepdown).first)
		{
			currentWeight += (*stepup).second;
			//cout << "UP " << (*stepup).second << "\t" << currentWeight << "\t" << (*stepup).first << "\n";
			stepup++;
			occupancy++;
		}
		else
		{
			currentWeight -= (*stepdown).second;
			//cout << "DOWN " << (*stepdown).second << "\t" << currentWeight << "\t" << (*stepdown).first << "\n";
			stepdown++;
			occupancy--;
		}
	}
#ifdef DBG_SPATIAL
	if (currentPos > cuts[curL+shift])
	{
		cout << "lost something " << previousPos << " " << currentPos << " " << cuts[curL+shift] << " " << curL << " " << shift << " " << rev <<  " " << n << "\n";
		for (int k=0 ; k<n+1; k++) {
			cerr << cuts[k] << "\n";
		}
	}
#endif
	length[curL] += cLength;
	Olength[curL] += cOLength;
	count[curL] += cCount;
	return;
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

	newIs->ISREF.occupancy = occupancy;

	// notify all tips on the current trajectory of insertion
	tip = notificationList.begin();
	while (tip != notificationList.end())
	{
		(**tip).notifyInsert(newIs);
		tip++;
	}


	return;
}

int Trajectory::numberOfSegmentsPassingPoint(double position)
{
  int occupancy = 0;

  list<Segment*>::iterator seg;
	seg = segments.begin();
	while (seg != segments.end())
	{
		(**seg).mt->updateLength();

		if ((**seg).crossesPoint(position))
		{	occupancy++;	}

		seg++;
	}
  return occupancy;
}
/*
int Trajectory::numberOfOtherSegmentsPassingIntersection(IntersectionItr& refis, Segment* refseg)
{
  int occupancy = 0;

  list<Segment*>::iterator seg;
	seg = segments.begin();
	while (seg != segments.end())
	{
    if (*seg != refseg)
    {
  		(**seg).mt->updateLength();
  
	  	if ((**seg).crossesIntersection(refis))
	  	{	occupancy++;	}
    }

		seg++;
	}
  return occupancy;
}
*/

// This version seems to be faster
int Trajectory::numberOfOtherSegmentsPassingIntersection(IntersectionItr& refis, Segment* refseg)
{
  int occupancy = 0;

	list<Segment*>::iterator segItr;
	IntersectionItr is;
  segItr = segments.begin();
  while (segItr != segments.end())
  {
    if ((*segItr) != refseg)
    {
      is = (**segItr).startItr;
      if ((**segItr).dir == ::forward)
      {
        while (++is != (**segItr).endItr)
        {
          if (is == refis)
          {
            occupancy++;
            break;
          }
        }
      }
      else
      {
        while (--is != (**segItr).endItr)
        {
          if (is == refis)
          {
            occupancy++;
            break;
          }
        }
      }
    }
    segItr++;
  }
  return occupancy;
}


int Trajectory::numberOfSegmentsDeflectingAtIntersection(IntersectionItr& is, Segment* refseg)
{
  int occupancy = 0;

  list<Segment*>::iterator seg;
	seg = segments.begin();
	while (seg != segments.end())
	{
    if (*seg != refseg)
    {
  		(**seg).mt->updateLength();
  
	  	if ((**seg).deflectsAtIntersection(is, refseg->dir))
	  	{	occupancy++;	}
    }

		seg++;
	}
  return occupancy;
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



/****************************** METATRAJECTORY FUNCTIONS **************************/

MetaTrajectory::MetaTrajectory(System* const s, double mDist, const SurfaceVector sv) : 
		system(s),
		startOfMetaTrajectory(sv),
		maxDistance(mDist),
		doesItIntersectOccupied(false)
{
	#ifdef DBG_ASTER
	cout << "DBG/ASTER: MetaTrajectory created\n";
	#endif

	addFirstTrajectory();
	#ifdef MTBASED_NUCLEATION_PROBABILITY
	crossedRegions.push_back(sv.region); 
	#endif

/*int regIDcounter;
for (regIDcounter=0;  regIDcounter<system->geometry->regions.size(); ++regIDcounter)
{
	if ( system->geometry->regions[regIDcounter]==sv.region ) 
	{	
		cout << "Aster is in region " << regIDcounter << endl;
		break;
	}
}*/



	bool noIntersection = true;
	
	while ( metaTrajectoryLength < maxDistance && noIntersection && !doesItIntersectOccupied )
	{
		
		TrajectoryVector tv = system->geometry->extendTrajectory(trajectoriesOfRay.back(), firstDirection);
/*int regIDcounter2;
for (regIDcounter2=0;  regIDcounter2<system->geometry->regions.size(); ++regIDcounter2)
{
	if ( system->geometry->regions[regIDcounter2]==tv.trajectory->base.region ) 
	{	
		cout << "First continuation is in region " << regIDcounter2 << endl;
		break;
	}
}
exit(1);*/
		trajectoriesOfRay.push_back(tv.trajectory);
		#ifdef MTBASED_NUCLEATION_PROBABILITY
		crossedRegions.push_back(tv.trajectory->base.region); 
		#endif
		int nOfIntersections = tv.trajectory->intersections.size();
		if ( nOfIntersections < 1 )
		{
			cerr << "Error in counting the number of intersections of a trajectory part of a metatrajectory, or I don't understsand how that thing works\n";
			exit(1);
		}
		
		if ( nOfIntersections > 1 )
			noIntersection = false;
			
		if ( noIntersection )
		{
			metaTrajectoryLength += tv.trajectory->length;
		}
		else 
		{
			IntersectionItr intItr = tv.trajectory->wallBegin();
			++intItr;
			while ( intItr != tv.trajectory->wallEnd() && !doesItIntersectOccupied )
			{
				if ( intItr->second.mirror->second.occupancy != 0/* occ != 0 */)                
				{
					nsc.numOfSegmentsInBundle = intItr->second.mirror->second.occupancy;                     //for debugging, will be removed
					double posOnTraj;
					if ( firstDirection==::forward )
					{
						posOnTraj = intItr->first;
					}
					else
					{
						posOnTraj = tv.trajectory->length - intItr->first;
					}
					metaTrajectoryLength += posOnTraj;               
					nsc.motherTrajectory = intItr->second.otherTrajectory;
					nsc.posOnMother = intItr->second.mirror->first;
					nsc.position = intItr->second.mirror->second.otherTrajectory->base;
					//nsc.position = intItr->second.otherTrajectory->base;
					nsc.position.region->translateVector( nsc.position, intItr->first );                            //here we get the SurfaceVector at the end of the meta trajectory (intersection with segment)
					doesItIntersectOccupied = true;
					noIntersection = false;
					break;
				}
				else
				{
					noIntersection = true;
					++intItr;
				}
			}
			
			if ( !doesItIntersectOccupied )
				metaTrajectoryLength += tv.trajectory->length;
			
		}
		
	} 

	return;	
}


MetaTrajectory::~MetaTrajectory()
{
	return;
}

void MetaTrajectory::addFirstTrajectory()
{
	TrajectoryVector tv = system->geometry->createTrajectory(startOfMetaTrajectory, false);
	firstDirection = tv.dir;
//	tv.dir=::forward;
//	firstDirection = tv.dir;
	trajectoriesOfRay.push_back(tv.trajectory);
	double angleMetaTrajectory = startOfMetaTrajectory.angle;//tv.trajectory->base.angle;//
//cout << "Angle first=" << angleMetaTrajectory*180/PI << ", with direction=" << tv.dir << endl;	
	SurfaceVector wallEndCoord;
	double wallX = (system->p.geomParam1)*0.5/(system->p.geomParam3);
	double wallY = (system->p.gridAspectratio)*(system->p.geomParam2)*0.5/(system->p.geomParam3);

	if ( angleMetaTrajectory == 0 )
	{
		wallEndCoord.x = wallX;
		wallEndCoord.y = startOfMetaTrajectory.y;
	}
	else if (  angleMetaTrajectory < PI/2 )                                    // TODO: correct for when angles are pi/2 or 3pi/2
	{
		if ( angleMetaTrajectory < atan( (wallY-startOfMetaTrajectory.y) / (wallX-startOfMetaTrajectory.x) ) )
		{
			wallEndCoord.x = wallX;
			wallEndCoord.y = tan( angleMetaTrajectory )*(wallX-startOfMetaTrajectory.x) + startOfMetaTrajectory.y;
		}
		else
		{
			wallEndCoord.x = ( wallY-startOfMetaTrajectory.y ) / tan( angleMetaTrajectory ) + startOfMetaTrajectory.x;
			wallEndCoord.y = wallY;
		}
	}
	else if ( angleMetaTrajectory == PI/2 )
	{
		wallEndCoord.x = startOfMetaTrajectory.x;
		wallEndCoord.y = wallY;
	}
	else if ( angleMetaTrajectory > PI/2 && angleMetaTrajectory < PI )
	{
		if ( angleMetaTrajectory-PI/2 < atan ( (wallX + startOfMetaTrajectory.x) / (wallY-startOfMetaTrajectory.y) ) )
		{
			wallEndCoord.x = -tan( angleMetaTrajectory-PI/2 )*(wallY-startOfMetaTrajectory.y) + startOfMetaTrajectory.x;
			wallEndCoord.y = wallY;
		}
		else
		{
			wallEndCoord.x = -wallX;
			wallEndCoord.y = (startOfMetaTrajectory.x + wallX) / tan( angleMetaTrajectory-PI/2 ) + startOfMetaTrajectory.y;
		}
	}
	else if ( angleMetaTrajectory == PI )
	{
		wallEndCoord.x = -wallX;
		wallEndCoord.y = startOfMetaTrajectory.y;
	}
	else if ( angleMetaTrajectory > PI && angleMetaTrajectory < 3*PI/2 )
	{
		if ( angleMetaTrajectory-PI < atan ( (wallY+startOfMetaTrajectory.y) / (wallX+startOfMetaTrajectory.x)) )
		{
			wallEndCoord.x = -wallX;
			wallEndCoord.y = tan( angleMetaTrajectory - PI )*(-wallX-startOfMetaTrajectory.x) + startOfMetaTrajectory.y;
		}
		else
		{
			wallEndCoord.x = ( -wallY-startOfMetaTrajectory.y ) / tan( angleMetaTrajectory - PI ) + startOfMetaTrajectory.x;
			wallEndCoord.y = -wallY;
		}
	}
	else if ( angleMetaTrajectory == 3*PI/2 )
	{
		wallEndCoord.x = startOfMetaTrajectory.x;
		wallEndCoord.y = -wallY;
	}
	else if ( angleMetaTrajectory > 3*PI/2 )
	{
		if ( angleMetaTrajectory-3*PI/2 < atan ( wallX-startOfMetaTrajectory.x / (wallY+startOfMetaTrajectory.y) ) )
		{
			wallEndCoord.x = ( wallY+startOfMetaTrajectory.y )*tan( angleMetaTrajectory - 3*PI/2 ) + startOfMetaTrajectory.x;
			wallEndCoord.y = -wallY;
		}
		else
		{
			wallEndCoord.x = wallX;
			wallEndCoord.y = -( wallX-startOfMetaTrajectory.x ) / tan( angleMetaTrajectory - 3*PI/2 ) + startOfMetaTrajectory.y;
		}
	}
	else
	{
		cout << angleMetaTrajectory*180/PI << endl;
		exit(42);
	}
	

	double firstMetaLength = sqrt( pow(wallEndCoord.x - startOfMetaTrajectory.x, 2) + pow(wallEndCoord.y - startOfMetaTrajectory.y, 2) );
	double posAsterOnTraj = tv.trajectory->length - firstMetaLength;
	
	int nOfIntersections = tv.trajectory->intersections.size();
	if ( nOfIntersections < 1 )
	{
		cerr << "Error in counting the number of intersections of a trajectory part of a metatrajectory, or I don't understsand how that thing works\n";
		exit(1);
	}
	if ( nOfIntersections != 1 )
	{
		IntersectionItr intItr = tv.trajectory->wallBegin();
		++intItr;
		while ( intItr != tv.trajectory->wallEnd() && !doesItIntersectOccupied )
		{
			double posOnTraj;
			if ( firstDirection==::forward )
			{
				posOnTraj = intItr->first;
			}
			else
			{
				posOnTraj = tv.trajectory->length - intItr->first;
			}
						
			if ( intItr->second.mirror->second.occupancy/*occ*/ != 0 &&  posOnTraj - posAsterOnTraj > ZERO_CUTOFF )                       
			{
				//cout << "Position of intersection on metatrajectory=" << posOnTraj << endl;
				//cout << "Position of aster on metatrajectory=" << posAsterOnTraj << endl;
				nsc.numOfSegmentsInBundle = intItr->second.mirror->second.occupancy;
				metaTrajectoryLength = posOnTraj - posAsterOnTraj;
				//cout << "Length of metaTrajectory=" << metaTrajectoryLength << endl;
				nsc.motherTrajectory = intItr->second.otherTrajectory;
				nsc.posOnMother = intItr->second.mirror->first;
				nsc.position = intItr->second.mirror->second.otherTrajectory->base;                                                      // <------ QUESTA E' CORRETTA
				nsc.position.region->translateVector( nsc.position, intItr->first );
				/*cout << "Centre of the aster=(" << startOfMetaTrajectory.x << "," << startOfMetaTrajectory.y << ")" << endl;
				cout << "First intersection at=(" << nsc.position.x << "," << nsc.position.y << ")" << endl;
				cout << "WallEnd coordinates=(" << wallEndCoord.x << "," << wallEndCoord.y << ")" << endl;
				cout << "!!!!!!!!!!!!!!!!!!!" << endl;
				if (angleMetaTrajectory*180/PI > -171) exit(1);*/
				doesItIntersectOccupied = true;
			}
			else
			{
				++intItr;
			}
		}
	}
	
	if ( !doesItIntersectOccupied )
		metaTrajectoryLength = firstMetaLength;
		
	return;

}


void MetaTrajectory::removeTrajectories() 
{
	vector <Trajectory*>:: iterator it = trajectoriesOfRay.begin();
	while ( it !=  trajectoriesOfRay.end() )
	{	
		(**it).conditionalRemove();
		 trajectoriesOfRay.erase(it);
	}
	return;
}
