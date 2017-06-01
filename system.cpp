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
#include <limits>

#include <float.h>		// necessary for DBL_MANT_DIG

#pragma warning(disable:981)
//#pragma warning(disable:177)
#pragma warning(disable:522)
#pragma warning(disable:383)

#include "MersenneTwister.h"
#include "corticalSim.h"

#ifdef NO_INLINE
#include "inline.cpp"
#endif

/***************************  SYSTEM FUNCTIONS ****************************/

System::System(char* parFile) : 
		EventDescriptorID(0), 
		eventID(0), 
		memUsage(0),
		totalLength(0),
	    seedsLeft(0),
		systemTime(0),
		systemTimeOffset(0),
		
		stopSignal(false),
		stopping(false),
		timeQueue(this, &System::identity, &System::identity),
		vPlusQueue(this, &System::vPlusToTime, &System::timeToVPlus),
		totalSEventCount(0),
		totalInvalidDEventCount(0),
		totalValidDEventCount(0),
		boundaryCrossingCount(0),
		totalZipperCount(0),
		totalCrossoverCount(0),
		totalInducedCatastropheCount(0),
		totalLengthSeveringCount(0),
		totalIntersectionSeveringCount(0),

		countSegments(0),
		countTrajectories(0),
		countIntersections(0),

		growingTipsSpecial(0),
		growingTipsNormal(0),

		dataSaved(false),
		histogramAverageCount(0),
		mtLengthHistogram(this),
		mtLifetimeHistogram(this),
		segAngleLengthHistogram(this),
		segAngleLifetimeHistogram(this),
		
		currentTimeTag(0),
		wallClockStartTime(time(0)),
		nextSnapshotEventTime(0),
		nextStatusEventTime(0),
		p(this),
		geometry(NULL)
/* NOTE: some compilers may complain about the presence of 'this' in the initializer
* list. This is no problem as long as the 'this' pointer is not used to access any of 
* the members during the construction of the object
*/
{
  // read parameters
  p.initialize(parFile);

  // perform initialization that is done only on creation
  switch(p.geometry)
  {
    case g_periodic:
      geometry = new Periodic(Coord2D(p.geomParam1, p.geomParam2), this);
      break;
    case g_grid:
      geometry = new Grid(Coord2D(p.geomParam1, p.geomParam2), static_cast<int>(p.geomParam3+0.4999), this);
      break;
    case g_wormhole:
      geometry = new Wormhole(Coord2D(p.geomParam1, p.geomParam2), this);
      break;
    case g_cylinder:
      geometry = new Cylinder(p.geomParam1, p.geomParam2, this);
      break;
    case g_gridcylinder:
      geometry = new GridCylinder(p.geomParam1, p.geomParam2, static_cast<int>(p.geomParam3+0.4999), this);
      break;
    case g_box:
      geometry = new Box(p.geomParam1, p.geomParam2, p.geomParam3, this);
      break;
    case g_pancake:
      geometry = new Pancake(p.geomParam1, this);
      break;
    default:
      cerr << "ERROR: no valid geometry handler present. aborting.\n";
      exit(-1);
  }

  randomGen.seed(static_cast<int>(p.seed));
  makeBinomialTable();


  seedsLeft += static_cast<int> (0.499 + p.preSeededSeedDensity * geometry->area) ;
	cout << "preSeeding "<< seedsLeft << " special nucleation events.\n";
	// Note that if forbiddenZones==1, a fraction of these events will be 'non-events'

  // perform further initialization
  mtLengthHistogram.setBins(p.hiresLengthHistogramBins);
  mtLifetimeHistogram.setBins(p.hiresLifetimeHistogramBins);
  segAngleLengthHistogram.setBins(p.loresAngleHistogramBins, p.loresLengthHistogramBins);
  segAngleLifetimeHistogram.setBins(p.loresAngleHistogramBins, p.loresLifetimeHistogramBins);


  initializeOutput();

  p.writeToFile();
  
  determineStochasticEvent();
  timeQueue.pushGlobal(p.measurementInterval, status);
  nextStatusEventTime = p.measurementInterval;
  if (p.movieEnabled)
  {
    timeQueue.pushGlobal(p.movieFrameInterval, snapshot);
    nextSnapshotEventTime = p.movieFrameInterval;
  }
  if (p.newParameterReadInterval > ZERO_CUTOFF)
  {
    timeQueue.pushGlobal(p.newParameterReadInterval, parameter_change);
    nextParameterEventTime = p.newParameterReadInterval;
  }

  return;
}


System::~System()
{


  //explicitly disassemble all microtubules
  growing_mts.removeAll();
  shrinking_mts.removeAll();


  return;	
}

void System::run(double simTime)
{
  int wallTimePollingCounter = 0;
  int memoryPollingCounter = 0;

  timeQueue.pushGlobal(simTime - systemTimeOffset, stop);
  while (!stopSignal)
  {
	  if (!stopping)
	  {
		if (totalLength > p.densityLimit*geometry->area)
		{
			//density limit reached; perform one final measurement, save a movie snapshot and send a stop signal	
			cerr << "Density limit (" << p.densityLimit << " /micrometer) reached. Stopping.\n";
			timeQueue.pushGlobal(systemTime, status);
			if (p.movieEnabled)
        			timeQueue.pushGlobal(systemTime, snapshot);
			timeQueue.pushGlobal(systemTime + 2*ZERO_CUTOFF, stop);
			stopping=true;
		}

		if(++wallTimePollingCounter == CLOCK_POLLING_INTERVAL)
		{
		  if (difftime(time(0),wallClockStartTime) >= p.wallClockLimit)
		  {
			//time's up; perform one final measurement, save a movie snapshot and send a stop signal	
			cerr << "Wall clock limit (" << p.wallClockLimit << " seconds) reached. Stopping.\n";
			timeQueue.pushGlobal(systemTime, status);
			if (p.movieEnabled)
        			timeQueue.pushGlobal(systemTime, snapshot);
			timeQueue.pushGlobal(systemTime + 2*ZERO_CUTOFF, stop);
			stopping=true;
		  }
		  wallTimePollingCounter = 0;
		}
		if (++memoryPollingCounter == MEMORY_POLLING_INTERVAL)
		{
		  if (estimateMemoryFootprint() > 1024*1024*p.memoryLimit)
		  {
			cerr << "Memory limit (approx. " << p.memoryLimit << " MB) reached. Stopping.\n";
			timeQueue.pushGlobal(systemTime, stop);
			stopping=true;
		  }
		  memoryPollingCounter = 0;
		}

		if (eventID > QUEUE_FLUSH_INTERVAL)
		{
	#ifdef DBG_QMGMT
		  cout << "Flushing queues to reduce numerical drift.\n";
	#endif
		  if(!integrityCheck()) //moved forward (used to be after flush and reload)
			exit(-1);
		  flushAndReload();
		  if(!integrityCheck())
			exit(-1);
	#ifdef DBG_QMGMT
		  cout << "Queue reloaded.\n";
	#endif
		}
	  }
    nextEvent();
  }

  writeMeasurementsToFile(0);
  closeFiles();
  return;
}

void System::emergencyBreak()
{
  writeMeasurementsToFile(0);	

  closeFiles();

	return;
}

int System::estimateMemoryFootprint()
{
	int total = 0;

	total += (growing_mts.size() + shrinking_mts.size())*(sizeof(Microtubule) + 2*sizeof(MTTip*));
	total += EventDescriptorMap.size() * (sizeof(pair<EventDescriptorIndex, EventDescriptor*>) + 3*sizeof(void*));
	total += countSegments * (sizeof(Segment) + sizeof(Segment*));
	total += countTrajectories * (sizeof(Trajectory));
	total += countIntersections * (2*(sizeof(pair<double, Intersection>) + 3*sizeof(void*)));

	total += (timeQueue.queue.size() + vPlusQueue.queue.size()) * (sizeof(DeterministicEvent) + 2*sizeof(void*));

	total += measurementHistory.size() * sizeof(Measurement);
	total += angleHistory.size() * p.angleHistogramBins * sizeof(double);
	total += angleNumberHistory.size() * p.discreteAngleNumber * sizeof(int);
	total += angleLengthHistory.size() * p.discreteAngleNumber * sizeof(double);

	return total;
}

bool System::integrityCheck()
{
	bool valid = true;
	
	if (!geometry->integrityCheck())
		valid = false;
	
	Microtubule* mt = growing_mts.first();
	while (mt != NULL)
	{
		if (!mt->integrityCheck())
			valid = false;
		mt = mt->next();
	}
	mt = shrinking_mts.first();
	while (mt != NULL)
	{
		if (!mt->integrityCheck())
			valid = false;
		mt = mt->next();
	}

	return valid;
}

void System::flushAndReload(bool refreshParameterEvent)
{
	double timeOffset = systemTime;

	// update all microtubule lengths.
	updateAll();
	// recalculate total lengths to avoid drift
	// in DBG mode, also check for large deviations
	double regionLength;
	double systemLength = 0;
	int ridx;
	Trajectory* tptr;

	for (ridx = 0; ridx < geometry->regions.size(); ridx++)
	{
		regionLength = 0;
		tptr = geometry->regions[ridx]->trajectories.first();
		while (tptr != NULL)
		{		
			regionLength += tptr->segmentLength();
			tptr = tptr->next();
		}
		if (abs(geometry->regions[ridx]->totalLength - regionLength) > ZERO_CUTOFF*max(100.,regionLength))
		{
			cerr << "Unacceptable drift in region length (drift = " << abs(geometry->regions[ridx]->totalLength - regionLength) << "). Exiting. [try lowering flush interval]\n";
			emergencyBreak();
			exit(-1);
		}
		geometry->regions[ridx]->totalLength = regionLength;
		systemLength += regionLength;
	}
	if (abs(totalLength - systemLength) > ZERO_CUTOFF*max(1.,totalLength))
	{
		cerr << "Unacceptable drift in system length (drift = " << abs(totalLength - systemLength) << "). Exiting. [try lowering flush interval]\n";
		emergencyBreak();
		exit(-1);
	}
	totalLength = systemLength;

	// reset master timer
	systemTimeOffset += timeOffset;
	systemTime = 0;
	
	timeQueue.flush();
	vPlusQueue.flush();
	EventDescriptorMap.clear();
	eventID=0;
	EventDescriptorID=0;

	// now, reload all events. This is a PROBLEM for global events that have not been saved. They need to be reinserted manually.
	determineStochasticEvent();
	timeQueue.pushGlobal(p.stopTime - systemTimeOffset, stop);
	nextStatusEventTime -= timeOffset;
	timeQueue.pushGlobal(nextStatusEventTime, status);

	if (refreshParameterEvent)
	{
		nextParameterEventTime -= timeOffset;
		timeQueue.pushGlobal(nextParameterEventTime, parameter_change);
	}

	if (p.movieEnabled)
	{
		nextSnapshotEventTime -= timeOffset;
		timeQueue.pushGlobal(nextSnapshotEventTime, snapshot);
	}
	
	Microtubule* mt;
	Segment* seg;
	mt = growing_mts.first();
	while (mt != NULL)
	{
		// this can go into a changetime function
		mt->nucleationTime -= timeOffset;
		seg = mt->segments.first();
		while (seg  != NULL)
		{
			seg->nucleationTime -= timeOffset;
			seg = seg->next();
		}

		mt->plus.event.index = registerEventDescriptor(&(mt->plus.event));
		mt->minus.event.index = registerEventDescriptor(&(mt->minus.event));
		mt->disappearEvent.index = registerEventDescriptor(&(mt->disappearEvent));

		mt->plus.determineEvent();
		mt->minus.determineEvent();
		mt->setDisappearEvent();
		
		mt = mt->next();
	}
	mt = shrinking_mts.first();
	while (mt != NULL)
	{
		mt->nucleationTime -= timeOffset;
		seg = mt->segments.first();
		{
			seg->nucleationTime -= timeOffset;
			seg = seg->next();
		}

		mt->plus.event.index = registerEventDescriptor(&(mt->plus.event));		// maybe put into member function
		mt->minus.event.index = registerEventDescriptor(&(mt->minus.event));		// maybe put into member function
		mt->disappearEvent.index = registerEventDescriptor(&(mt->disappearEvent));		// maybe put into member function

		mt->plus.determineEvent();
		mt->minus.determineEvent();
		mt->setDisappearEvent();

		mt = mt->next();
	}

	
	return;
}

#ifndef NO_INLINE 
inline 
#endif 
double LambertW(double x)
{
/*
Function is based on the numerical procedure described in 

Real Values of the W-Function, Barry et al., ACM Trans. on Math. Software, Vol. 21, No. 2, June 1995, pp. 161--171.

Only a single convergence pass is used, yielding 16 significant digits.

*/

	double out;

	const double expMin1 = exp(-1.0);
	const double zeroCutoff = pow(2.0, -DBL_MANT_DIG/6.0);
	const double etaCutoff = 34.0 * pow(2.0, 2.0*(1.0-DBL_MANT_DIG)/7.0);

	if (x < -expMin1) return -VERY_LARGE;
	
	if (abs(x) < zeroCutoff)
	{
		out = x*(60 + 114*x + 17*x*x)/(60 + 174*x + 101*x*x);
	} 
	else if (x > 20.0)
	{
		double h = exp(- 1.124491989777808 /(0.4225028202459761 + log(x)));
		out = log(x /log(x/pow(log(x),h)));

		double z_n = log(x/out) - out;
		double temp = 2.0*(1.0 + out)*(1.0 + out + z_n*2.0/3.0);
		double e_n = (z_n/(1.0 + out))*(temp - z_n)/(temp - 2.0*z_n);
		out = out*(1.0 + e_n);
	}
	else
	{

		double eta = 2.0 + 2.0*exp(1.0)*x;
		double sqrt_eta = sqrt(eta);
		double D;
		if (sqrt_eta >= etaCutoff)
		{
			double n2 = 4.612634277343749*sqrt(sqrt(sqrt_eta + 1.09556884765625));
			double n1 = (4.0 - 3.0*sqrt(2.0) + (2.0*sqrt(2.0) - 3.0)*n2)/(sqrt(2.0) - 2.0);
			D = n1*sqrt_eta/(n2 + sqrt_eta);

			out = -1.0 + sqrt_eta/(1.0 + sqrt_eta/(3.0 + D));

			double z_n = log(x/out) - out;
			double temp = 2.0*(1.0 + out)*(1.0 + out + z_n*2.0/3.0);
			double e_n = (z_n/(1.0 + out))*(temp - z_n)/(temp - 2.0*z_n);
			out = out*(1.0 + e_n);
		}
		else
		{
			D = sqrt_eta/(8.0/3.0 + sqrt_eta/(135.0/83.0 + sqrt_eta/(166.0/39.0 + sqrt_eta*3167.0/3549.0)));
			out = -1.0 + sqrt_eta/(1.0 + sqrt_eta/(3.0 + D));
		}
	}

#ifdef DBG_ACID_TEST
	double test = (x - out*exp(out))/((1.0 + out)*exp(out));
	if (abs(test) > ZERO_CUTOFF*abs(x))
		cerr << "LambertW() accuracy [" << test << "] below threshold. Output may contain inaccuracies.\n";
#endif

	return out;
}


double System::vPlusToTime(double dist)
{
    return dist;

}

double System::timeToVPlus(double time)
{
    return time;

}


void System::advanceTime(double newTime)
{
	#ifdef DBG_SYSTEM
	cout << "DBG/SYSTEM: advancing time by " << newTime - systemTime << " seconds. New system time is " << newTime << ".\n";
	#endif

	timeQueue.advanceTime(newTime);
	vPlusQueue.advanceTime(newTime);
	systemTime = newTime;

	totalLength += timeQueue.progression(currentTimeTag)*(p.vMin*shrinking_mts.size() - p.vTM*(growing_mts.size() + shrinking_mts.size()) ) 
				 + vPlusQueue.progression(currentTimeTag)*(p.vPlus*growing_mts.size());

	if (++currentTimeTag == POSITION_CACHE_SIZE)
	{
		currentTimeTag = 0;
		updateAll(true);	// updates lengths and forces all tags to tag=0.
	}

	// Note: the following statements should be performed *after* updateAll(true)!
	timeQueue.storeTime(currentTimeTag);
	vPlusQueue.storeTime(currentTimeTag);

	return;
}


void System::updateAll(bool forceUpdate)
/* Updates the lengths of all regions and microtubules. Convenient for snapshots, etc. 
*/
{
	for (int ridx = 0; ridx < geometry->regions.size(); ridx++)
		geometry->regions[ridx]->updateRegionLength(forceUpdate);

	Microtubule* mt;
	mt = growing_mts.first();
	while (mt!= NULL)
	{
		mt->updateLength(forceUpdate);
		
		#ifdef DBG_ACID_TEST
		if (mt->segments.size() == 1)
		{
			if ((mt->segments.first()->end - mt->segments.first()->start)*mt->segments.first()->dir < -ZERO_CUTOFF)
			{
				cerr << "ERROR: negative segment length after updating.\n";
				exit(-1);
			}
		}
		#endif
		
		mt = mt->next();
	}
	mt = shrinking_mts.first();
	while (mt!= NULL)
	{
		mt->updateLength(forceUpdate);

		#ifdef DBG_ACID_TEST
		if (mt->segments.size() == 1)
		{
			if ((mt->segments.first()->end - mt->segments.first()->start)*mt->segments.first()->dir < -ZERO_CUTOFF)
			{
				cerr << "ERROR: negative segment length after updating.\n";
				exit(-1);
			}
		}
		#endif
		mt = mt->next();
	}
	
	return;
}
	
void System::nextEvent()
{


	double tqValue = numeric_limits<double>::max();
	double dqValue = numeric_limits<double>::max();

	// maybe clean up a little
	if (!timeQueue.empty())
		tqValue = timeQueue.firstEventTime();
	if (!vPlusQueue.empty())
		dqValue = vPlusQueue.firstEventTime();
	
	if ((nextStochasticEventTime < tqValue) && (nextStochasticEventTime < dqValue))
	{
		// execute stochastic event

		totalSEventCount++;
		#ifdef DBG_EVENT
		cout << "DBG/Event: executing stochastic event [time=" << nextStochasticEventTime << ", type=" << nextStochasticEventType << "]\n";
		#endif

		advanceTime(nextStochasticEventTime);
		
		// execute next stochastic event
		switch(nextStochasticEventType)
		{
		case catastrophe:
			handleCatastropheEvent();
			break;
		case rescue:
			handleRescueEvent();
			break;
		case katanin:
			handleSeveringEvent();
			break;
		case severingAtCross:
			handleSeveringAtCrossEvent();
			break;
		case nucleation:
			handleNucleationEvent();
			break;
		case preSeededNucleation:
			handleNucleationEvent(true);
			break;
		}	
		determineStochasticEvent();
	}
	else
	{
		double nextEventTime;
		// deterministic event
		DeterministicEvent event;
		if (tqValue < dqValue)
		{
			event = timeQueue.pop();
			nextEventTime = tqValue;
		}
		else
		{
			event = vPlusQueue.pop();
			nextEventTime = dqValue;
		}

		if (event.infoIdx != -1)
		{
			map<EventDescriptorIndex, EventDescriptor*>::iterator eventItr = EventDescriptorMap.find(event.infoIdx);
			// Invalidated events are not removed from the queue, but the nextEvent for the segment *is* updated.
			// Therefore, check whether the tip still exists and whether the tag still matches
			if ((eventItr != EventDescriptorMap.end()) 
				&& (eventItr->second->tag == event.tag)) 
			{
				#ifdef DBG_EVENT
				cout << "DBG/Event: executing deterministic event [time=" << event.eventTimeDist << ", type=" << eventItr->second->type << ", tag=" << event.tag << "]\n";
				#endif

				totalValidDEventCount++;
				advanceTime(nextEventTime);
				eventItr->second->mt->handleEvent(eventItr->second);
				determineStochasticEvent();
			}
			else // no longer a valid event
			{
				totalInvalidDEventCount++;
				#ifdef DBG_EVENT
				cout << "DBG/Event: Event is no longer valid. [tag= " << event.tag << "]\n";
				#endif
			}
		}
		else
		{
			#ifdef DBG_EVENT
			cout << "DBG/Event: executing global event [time=" << event.eventTimeDist << ", type=" << event.global_type << "]\n";
			#endif
	
			// execute next global event
			advanceTime(event.eventTimeDist);

			handleGlobalEvent(event);
		}
	}

	
	return;
}

#ifndef NO_INLINE
inline 
#endif
void System::handleCatastropheEvent()
{
	
		int tipNumber, totalRelevantTips;
		int ridx = 0;
		RegionParameterType paramType;
		RegionMTTipTag tipTag;

		// determine whether we're looking for a 'normal' or 'special' tip
		if (randomGen.rand() < (double)growingTipsNormal/((double)growingTipsNormal + (double)p.catastropheMultiplier*growingTipsSpecial))
		{
			// catastrophe takes place on a regular face of the geometry (kCat)
			paramType = r_param_normal;
			totalRelevantTips = growingTipsNormal;
		}
		else
		{
			// catastrophe takes place on a special face of the geometry (kCat*catastropheMultiplier)
			paramType = r_param_modified;
			totalRelevantTips = growingTipsSpecial;
		}

		// select the number of the tip that will undergo a catastrophe
		tipNumber = randomGen.randInt(totalRelevantTips - 1);

		// determine the region this tip is located in
		// depending on whether it's a high or low number, start scanning from the beginning or end
		if (tipNumber < totalRelevantTips/2)
		{
			while (true)
			{
				if (geometry->regions[ridx]->parameterType == paramType)
				{
					if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size())
						break;
					tipNumber -= geometry->regions[ridx]->growingPlusTipList.size();
				}
				ridx++;
			}
		}
		else
		{
			// flip the direction of indexing: back to front
			tipNumber = totalRelevantTips- 1 - tipNumber;
			ridx = geometry->regions.size()-1;
			while (true)
			{
				if (geometry->regions[ridx]->parameterType == paramType)
				{
					if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size())
						break;
					tipNumber -= geometry->regions[ridx]->growingPlusTipList.size();
				}
				ridx--;
			}
			// flip the direction of the result for the next step
			tipNumber = geometry->regions[ridx]->growingPlusTipList.size() - 1 - tipNumber;
		}

		// Now, select the tip with # 'tipNumber' from the current region.
		// Again, if the number is high within the number of tips in the region, start from the end
		if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size()/2)
		{
			tipTag = geometry->regions[ridx]->growingPlusTipList.begin();
			while (tipNumber >0)
			{
				tipTag++;
				tipNumber--;
			}
		}
		else
		{
			tipTag = --(geometry->regions[ridx]->growingPlusTipList.end());
			while (tipNumber < geometry->regions[ridx]->growingPlusTipList.size()-1)
			{
				tipTag--;
				tipNumber++;
			}
		}



		(*tipTag)->mt->catastrophe();

		return;
}



#ifndef NO_INLINE
inline 
#endif
void System::handleRescueEvent()
{
		int ridx = 0;
		int tipNumber;
		RegionMTTipTag tipTag;

		// select the number of the tip that will be rescued
		tipNumber = randomGen.randInt(shrinking_mts.size() - 1);

		// determine the region this tip is located in
		// depending on whether it's a high or low number, start scanning from the beginning or end
		if (tipNumber < shrinking_mts.size()/2)
		{
			while (tipNumber >= geometry->regions[ridx]->shrinkingPlusTipList.size())
			{
				tipNumber -= geometry->regions[ridx]->shrinkingPlusTipList.size();
				ridx++;
			}
		}
		else
		{
			// flip the direction of indexing: back to front
			tipNumber = shrinking_mts.size()- 1 - tipNumber;
			ridx = geometry->regions.size()-1;
			while (tipNumber >= geometry->regions[ridx]->shrinkingPlusTipList.size())
			{
				tipNumber -= geometry->regions[ridx]->shrinkingPlusTipList.size();
				ridx--;
			}
			// flip the direction of the result for the next step
			tipNumber = geometry->regions[ridx]->shrinkingPlusTipList.size() - 1 - tipNumber;
		}

		// Now, select the tip with # 'tipNumber' from the current region.
		// Again, if the number is high within the number of tips in the region, start from the end
		if (tipNumber < geometry->regions[ridx]->shrinkingPlusTipList.size()/2)
		{
			tipTag = geometry->regions[ridx]->shrinkingPlusTipList.begin();
			while (tipNumber >0)
			{
				tipTag++;
				tipNumber--;
			}
		}
		else
		{
			tipTag = --(geometry->regions[ridx]->shrinkingPlusTipList.end());
			while (tipNumber < geometry->regions[ridx]->shrinkingPlusTipList.size()-1)
			{
				tipTag--;
				tipNumber++;
			}
		}


		// call the rescue function for the selected tip
		(*tipTag)->mt->rescue();

		return;


}


#ifndef NO_INLINE
inline 
#endif
void System::handleSeveringEvent()
{
	double cutLength;

	Segment* cutSeg;

	randomPositionOnMicrotubule(cutLength, cutSeg);
	if (cutSeg != NULL)	
		cutSeg->mt->sever(cutSeg,cutLength);
	else
		cerr << "ERROR (non-fatal): severing position outside total MT length. Ignoring event.\n";


	return;
}

void System::randomPositionOnMicrotubule(double& randomPos, Segment*& randomSeg)
{
	int ridx = 0;

	// select position at which severing will take place
	double cutLength = randomGen.randDblExc(totalLength);


	// first, select relevant region. Go from both sides to reduce time.
	if (cutLength < 0.5*totalLength)
	{
		ridx = 0;
		while (true)
		{
			geometry->regions[ridx]->updateRegionLength();
			if (cutLength < geometry->regions[ridx]->totalLength)
				break;
			cutLength -= geometry->regions[ridx]->totalLength;
			ridx++;
			if (ridx == geometry->regions.size())
			{
				// apparently, cutLength was *just* too long...
				cutLength -= ZERO_CUTOFF;
#ifdef DBG_ACID_TEST
				if (cutLength > geometry->regions[ridx-1]->totalLength)
				{
					cerr << "ERROR: cannot locate intersection location\n";
					exit(-1);
				}
#endif
				break;
			}
		}
	}
	else
	{
		ridx = geometry->regions.size()-1;
		cutLength = totalLength - cutLength;
		while (true)
		{
			geometry->regions[ridx]->updateRegionLength();
			if (cutLength < geometry->regions[ridx]->totalLength)
				break;
			cutLength -= geometry->regions[ridx]->totalLength;
			ridx--;
			if (ridx == -1)
			{
				// apparently, cutLength was *just* too long...
				cutLength -= ZERO_CUTOFF;
#ifdef DBG_ACID_TEST
				if (cutLength > geometry->regions[0]->totalLength)
				{
					cerr << "ERROR: cannot locate intersection location\n";
					exit(-1);
				}
#endif
				break;
			}
		}
		//and wrap it back
		cutLength = geometry->regions[ridx]->totalLength - cutLength;
	}

	// Now that the area has been selected, cycle over trajectories and their segments

	Trajectory* trptr;
	list<Segment*>::iterator seg;
	double temp;
	bool stopFlag = false;

	// Again, let the search direction depend on the expected position
	if (cutLength < 0.5*(geometry->regions[ridx]->totalLength))
	{
		// search in the forward direction. Start with the first trajectory
		trptr = geometry->regions[ridx]->trajectories.first();
		while (true)
		{
			// get the first segment on this trajectory
			seg = trptr->segments.begin();
			// and loop over all segments
			while (seg != trptr->segments.end())
			{
				// make sure the lengths are accurate
				(*seg)->mt->updateLength();
				// break when the correct segment is found
				if (cutLength < (temp = (*seg)->length()))
				{
					stopFlag = true;
					break;
				}
				else
				{
					// if not, continue
					cutLength -= temp;
					seg++;
				}
			}

			if (stopFlag)
				break;

			// if we ran out of trajectories, somehing is wrong. This could happen due to numerical inaccuracies.
			if (trptr->next() == NULL)
			{
				cerr << "Rare event: cutting on the edge\n";
				seg--;
				cutLength -= ZERO_CUTOFF;
				break;
			}
			trptr = trptr->next();
		}
	}
	else
	{
		// flip the length and trajectory search orders
		cutLength = geometry->regions[ridx]->totalLength - cutLength;
		trptr = geometry->regions[ridx]->trajectories.last();
		while (true)
		{
			// also, start with the last segment (one step beyond, in this case)
			seg = (trptr->segments.end());
			do
			{
				seg--;
				// make sure that the segment length is correct
				(*seg)->mt->updateLength();
				// break if the requested position lies on this segment
				if (cutLength < (temp = (**seg).length()))
				{
					stopFlag = true;
					break;
				}
				else
					cutLength -= temp;
			}
			while (seg != trptr->segments.begin());

			if (stopFlag)
				break;

			if (trptr->previous() == NULL)
			{
				cerr << "Rare event: cutting on the edge\n";
				seg++;
				cutLength -= ZERO_CUTOFF;
				break;
			}
			trptr = trptr->previous();
		}
		cutLength = (*seg)->length() - cutLength;
	}

	randomSeg = *seg;
	randomPos = cutLength;
	return;
}

#ifndef NO_INLINE
inline 
#endif
// TODO: verify this function
void System::handleSeveringAtCrossEvent()
	{
		Microtubule* cutMT;
		IntersectionItr cutIS;
		OccupiedIntersection* cutOccIS;
		Segment* cutSeg;
		int cutSegNumber;
		list<Segment*>::iterator testSeg;
		list<Segment*>::iterator finalTestSeg; 

		cutSeg = NULL; 
		cutOccIS = OccupiedIntersectionList.ElementAddress(randomGen.randInt(OccupiedIntersectionList.size()-1));
		if (p.crossSeveringTop || (randomGen.randInt(1)==0))
			// if MT (bundle) on top should be cut, select that one, or (if random) pick it with 50% chance.
			// Note that no effort is made to ensure proportional selection for bundles of different thickness.
			cutIS = cutOccIS->intersectionToCut; 
		else
			cutIS = cutOccIS->intersectionToCut->second.mirror;
		double diffAngle = cutIS->second.otherTrajectory->base.region->intersectionAngle(
				cutIS->second.otherTrajectory,
				cutIS->second.mirror->second.otherTrajectory);
		if (diffAngle > 0.5*PI)
			diffAngle = PI - diffAngle;
		// don't do anything if the angle between MTs is smaller than the cutoff angle
		if (diffAngle <= (p.crossSeveringStartAngle)*PI/180.)
			return;

		/*
		IMPORTANT NOTE: DO NOT REMOVE OCCUPIED INTERSECTION HERE! IT WILL BE REMOVED AUTOMATICALLY IN A LATER EVENT 
		cutIS->second.occupancy--;		// do not decrease: will happen very soon in a new event
		{
			removeOccupiedIntersection(cutIS->second);
		}	
		*/
		// find the segment to be cut
		cutSegNumber = randomGen.randInt(cutIS->second.occupancy - 1)  ; // nummers 0 - occupancy -1 . 
		if (cutSegNumber < cutIS->second.occupancy/2)	// search from the beginning. 
		// In principle, all searches could start from the beginning. Also starting from the end should save (on average) half of the time 
		{
			testSeg = cutIS->second.mirror->second.otherTrajectory->segments.begin();
			finalTestSeg = cutIS->second.mirror->second.otherTrajectory->segments.end();
			cutSegNumber++; 
			while (testSeg != finalTestSeg) // NOTE: test should not be necessary "while(true)"
			{
				if ((**testSeg).crossesIntersection(cutIS))
				{				
					cutSegNumber--;
					if (cutSegNumber == 0)
					{
						cutSeg = *testSeg;
						break;
					}
				}
				testSeg++;
			}
		}
		else	// search from the end
		{
			testSeg = cutIS->second.mirror->second.otherTrajectory->segments.end();
			finalTestSeg = cutIS->second.mirror->second.otherTrajectory->segments.begin();
			while (testSeg != finalTestSeg)
			{
				testSeg--;
				if ((**testSeg).crossesIntersection(cutIS))
				{				
					cutSegNumber++;
					if (cutSegNumber == cutIS->second.occupancy)
					{
						cutSeg = *testSeg;
						break;
					}
				}
			}
		}
		if(cutSeg != NULL)
		{
			cutMT = cutSeg->mt;
			cutMT->updateLength();
			if (cutMT != NULL)	
				cutMT->severAtCross(cutIS, cutSeg);
			else
				cerr << "ERROR (non-fatal): severing position outside total MT length. Ignoring event.\n";
		}
		else
			cerr << "ERROR (non-fatal): Did not find the right segment. Ignoring event.\n";


		return;
}

#ifndef NO_INLINE
inline 
#endif
void System::handleNucleationEvent(bool preSeeded)
{

	TrajectoryVector tv;
	SurfaceVector sv = geometry->randomSurfaceVector();
	Segment* nucSeg;
  NucleationType thisNuc;
	double pos,posOnSeg;
	double overrideAngle, rand,rand2, preTheta, cosTheta, temp;

	if (preSeeded)
	{
		--seedsLeft;
	    thisNuc = p.preSeededType;
	}
	else
	{
		thisNuc = p.nucleationType;
	}

	if ((p.forbiddenZones) && (sv.region->parameterType == r_param_modified))
		return;

  switch (thisNuc)
	{
		case nuc_isotropic:
			overrideAngle = sv.angle;
			break;
		case nuc_discreteAngles:
    
			overrideAngle = p.nucleationAngles[randomGen.randInt(p.discreteAngleNumber-1)] + PI*randomGen.randInt(1);
			break;
		case nuc_ellipse:
    		if (randomGen.randDblExc() >= p.nucleationHalfIsotropicDensity/(totalLength/geometry->area + p.nucleationHalfIsotropicDensity))
			{
				randomPositionOnMicrotubule(posOnSeg, nucSeg);
				pos = nucSeg->start + nucSeg->dir*posOnSeg;

				rand = randomGen.randDblExc();
				if (p.ellipseForwardAlongMT && rand >= p.ellipseLeftFraction + p.ellipseRightFraction)
				{
					tv.pos=pos;
					tv.trajectory = nucSeg->trajectory;
					if (rand >= 1. - p.ellipseBackwardFraction)
					{
						// Backward nucleation
						tv.dir = (nucSeg->dir == ::forward) ? backward : ::forward ;
					}
					else
					{
						// Forward nucleation
						tv.dir = nucSeg->dir;
					}

				}
				else
				{
					sv = nucSeg->trajectory->base;
					sv.region->translateVector(sv, pos);
					preTheta = randomGen.randDblExc()*2.*PI;
					cosTheta = cos(preTheta);
					rand2 = randomGen.randDblExc();
					rand2 = sqrt(rand2);
					temp = rand2*sqrt((1.-p.ellipseEpsilon*p.ellipseEpsilon)*(1.-cosTheta*cosTheta));
					if (preTheta > PI)
					{
						temp = - temp;
					} 
					overrideAngle = atan2(temp, (p.ellipseEpsilon + rand2*cosTheta));

					// rotate by +/- 40 degrees if wanted			
					preTheta = 0.;
					if ( rand < p.ellipseLeftFraction + p.ellipseRightFraction ) 
					{
						preTheta = PI/180.*p.ellipseSidewaysAngle;
						if (rand >= p.ellipseLeftFraction )
						{
							preTheta = -preTheta;
						}
					}
					else if ( rand >= 1.- p.ellipseBackwardFraction )
					{
						preTheta = PI;
					}
					overrideAngle += preTheta;
					sv.angle = nucSeg->trajectory->base.angle + overrideAngle;
					if (nucSeg->dir == backward)
						sv.angle += PI;

					// for safety: bring back to interval [0..2 PI]
					while (sv.angle > 2*PI)
						sv.angle -= 2*PI;
					while (sv.angle < 0)
						sv.angle += 2*PI;
#ifdef DBG_ASSERT
					if (sv.angle > 2*PI)
						cerr << "ERROR: function System::handleNucleationEvent: angle > 2*PI (" << sv.angle << ")\n";
					if (sv.angle < 0.)
						cerr << "ERROR: function System::handleNucleationEvent: angle < 0 (" << sv.angle << ")\n";
#endif


					sv.region->translateVector(sv, ZERO_CUTOFF);
					tv = geometry->createTrajectory(sv);
				}
				if ( true ) // extreem smerige constructie om gezeik van de compiler te voorkomen :-((
				{
					Microtubule* mt = growing_mts.create(this, tv);			
				}
				return; 
			}
			else
			{
				if ( p.ellipseReducedFreeRate == 0 || randomGen.randDblExc() < p.ellipseReducedFreeRateAcceptFraction ) {
					if ( abs(p.nucleationAlpha) < ZERO_CUTOFF ) 
					{
						overrideAngle = sv.angle;
						break;
					}
					// Lacking any break statement on purpose: control should flow to nuc_biased for p.nucleationAlpha > ZERO_CUTOFF
				}	
				else
					// nucleation event rejected.
					return;
			}
			// WARNING : do not insert statements here!
		case nuc_biased:		// preferentially at 1/4 pi + k*1/2 pi 
			// WARNING : nuc_biased should directly follow nuc_ellipse
			static int firstcall =  1;
			static double eps, fa1, fa2, gmax, binsize, invbinsize;
			static double *chi, *gam0;	
			static int firstcall_ps =  1;
			static double eps_ps, fa1_ps, fa2_ps, gmax_ps, binsize_ps, invbinsize_ps;
			static double *chi_ps, *gam0_ps;	
			double alpha;

			if (sv.region->accountingType == r_accounting_dontcount)
			{
				overrideAngle = sv.angle;
				break;
			}
      
  
			int bins = NUCLEATION_DISCRETIZATION_STEPS;
			int i;
			double gam, fa;


			if (preSeeded)
			{
				if (firstcall_ps)
				{
					if (abs(p.preSeededAlpha) < ZERO_CUTOFF) 
						cerr << "preSeededAlpha too small - this causes segmentation fault!\n";
					eps_ps = exp(- p.preSeededAlpha)/(exp(p.preSeededAlpha) - exp(- p.preSeededAlpha));
					fa1_ps = exp(p.preSeededAlpha);
					fa2_ps = exp(p.preSeededAlpha) - exp(- p.preSeededAlpha);
					chi_ps = new double[bins+1];
					gam0_ps = new double[bins+1];
					gmax_ps = 4./3. + 2*eps_ps;
					binsize_ps = gmax_ps/bins;
					invbinsize_ps = bins/gmax_ps;
					for (i=0; i<=bins; i++)		// table contains values for both 0 and gmax  
					{
						gam = binsize_ps * i; 
						gam0_ps[i] = gam;
						chi_ps[i] = 2*sqrt(1+eps_ps)*cos(1./3.*(2*PI-acos((1+1.5*eps_ps-1.5*gam)/pow((1+eps_ps),1.5))));		
					}
					firstcall_ps = 0;
				}
			}
			else 
			{
				if (firstcall)
				{
					if (abs(p.nucleationAlpha) < ZERO_CUTOFF) 
						cerr << "nucleationAlpha too small - this causes segmentation fault!\n";
					eps = exp(- p.nucleationAlpha)/(exp(p.nucleationAlpha) - exp(- p.nucleationAlpha));
					fa1 = exp(p.nucleationAlpha);
					fa2 = exp(p.nucleationAlpha) - exp(- p.nucleationAlpha);
					chi = new double[bins+1];
					gam0 = new double[bins+1];
					gmax = 4./3. + 2*eps;
					binsize = gmax/bins;
					invbinsize = bins/gmax;
					for (i=0; i<=bins; i++)		// table contains values for both 0 and gmax  
					{
						gam = binsize * i; 
						gam0[i] = gam;
						chi[i] = 2*sqrt(1+eps)*cos(1./3.*(2*PI-acos((1+1.5*eps-1.5*gam)/pow((1+eps),1.5))));		
					}
					firstcall = 0;
				}
			}

			while (1)
			{
				if (preSeeded) 
				{
					gam = randomGen.randDblExc(gmax_ps);
					i = static_cast<int>(gam/binsize_ps);
					gam = invbinsize_ps*((gam-gam0_ps[i])*chi_ps[i+1] + (gam0_ps[i+1]-gam)*chi_ps[i]);
					fa = fa1_ps - fa2_ps*gam*gam;
					alpha = p.preSeededAlpha;
				}
				else 
				{
					gam = randomGen.randDblExc(gmax);
					i = static_cast<int>(gam/binsize);
					gam = invbinsize*((gam-gam0[i])*chi[i+1] + (gam0[i+1]-gam)*chi[i]);
					fa = fa1 - fa2*gam*gam;
					alpha = p.nucleationAlpha;
				}
				if (randomGen.randDblExc(fa) < exp(alpha*cos(PI*gam)))
				{
					overrideAngle = 0.25*PI*(1+gam)+0.5*PI*randomGen.randInt(3);
					break;
				}
			}
			break;

	}	

	sv.angle = overrideAngle;

	tv = geometry->createTrajectory(sv);
	Microtubule* mt = growing_mts.create(this, tv);

	return;
}




void System::determineStochasticEvent()
{

	double nextEventInterval;
	double invTotalRate;
	double typeSelector;
	double logRandom;
	double temp;

    double constRate;
    double slopeRate;

    constRate = seedsLeft*p.preSeededRate +p.kNuc*geometry->area + p.kRes*shrinking_mts.size()\
            + p.kSev*totalLength
            + p.kCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial)
            + p.kCross*OccupiedIntersectionList.size();
    slopeRate = p.kSev*((p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size());

    logRandom = -log(randomGen.randDblExc()); //random number larger than 0 with an exponential distribution
    if (abs(slopeRate) < ZERO_CUTOFF)
    {
        nextEventInterval = logRandom / constRate;
    }
    else
    {
        if ((temp = pow(constRate,2) + 2*slopeRate*logRandom) < 0)
            nextEventInterval = VERY_LARGE; // select a time that is always larger than the disappearance time of all MTs
        else
            nextEventInterval = (-constRate + sqrt(temp))/slopeRate;
    }


    invTotalRate = 1.0/(constRate + slopeRate*nextEventInterval); //NOTE: will be negative for the case without stochastic events

#ifdef DBG_ACID_TEST
	if ((nextEventInterval < 0) || (nextEventInterval != nextEventInterval ))
	{
		cerr << "Invalid deterministic event time calculated. Exiting\n";
		exit(-1);
	}
#endif
	nextStochasticEventTime = systemTime + nextEventInterval;

	// note that in the case of a non-event (interval = VERY_LARGE), invTotalRate can be negative.
	// This does not matter, as long as there is a fallback event type selected (in this case katanin). Due to the timing, the
	// event will never be executed.
	typeSelector = randomGen();
	if (typeSelector < p.kNuc*geometry->area*invTotalRate)
		nextStochasticEventType = nucleation;
	else if (typeSelector < (p.kNuc*geometry->area + p.kCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial))*invTotalRate)
		nextStochasticEventType = catastrophe;
	else if (typeSelector < (p.kNuc*geometry->area + p.kCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + p.kRes*shrinking_mts.size())*invTotalRate)
		nextStochasticEventType = rescue;
	else if (typeSelector < (p.kNuc*geometry->area + p.kCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + p.kRes*shrinking_mts.size() + seedsLeft*p.preSeededRate)*invTotalRate)
		nextStochasticEventType = preSeededNucleation;
	else if (typeSelector < (p.kNuc*geometry->area + p.kCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + p.kRes*shrinking_mts.size() + seedsLeft*p.preSeededRate + p.kCross*OccupiedIntersectionList.size())*invTotalRate)
		nextStochasticEventType = severingAtCross;
	else
		nextStochasticEventType = katanin;

#ifdef DBG_EVENT
	cerr << "determined next stochastic event. type=" << nextStochasticEventType << ", time=" << nextStochasticEventTime << "\n";
#endif

	return;
}


void System::handleGlobalEvent(DeterministicEvent& event)
{

	updateAll();

	switch (event.global_type)
	{
		case status:
#ifndef NO_INTEGRITY_CHECK
			if (!integrityCheck())
			{
				cerr << "Integrity check failed. Exiting.\n";
				exit(-1);
			}
#endif

			performMeasurement();
			nextStatusEventTime = systemTime + p.measurementInterval;
			timeQueue.pushGlobal(nextStatusEventTime, status);

			break;	
		case snapshot:
			movieFile << "time " << systemTime + systemTimeOffset << "\n";
			OrderParameters order;
			geometry->getOrderParameters(order);
			movieFile << "order " << order.R << " " << order.Rdirector[0] << " " << order.Rdirector[1] << " " << order.Rdirector[2] << "\n";
			geometry->outputSnapshot(movieFile);
			movieFile << "endOfFrame\n";
			nextSnapshotEventTime = systemTime + p.movieFrameInterval;
			timeQueue.pushGlobal(nextSnapshotEventTime, snapshot);

			break;
		case stop:
			stopSignal = true;
			break;
		case parameter_change:
			cout << "Changing parameters. Reading from file [" << p.newParameterFile << "]\n";

			// reload

			flushAndReload(false);		// check out what to do with parameter change events.

			string oldDir = p.outputDir;
			if (p.reinitialize(p.newParameterFile.c_str()))
			{
				if (p.outputDir != oldDir)
				{
					cout << "Switching directories.\n";
					writeMeasurementsToFile(0);		// flush all measurements
					closeFiles();					// close all files
					initializeOutput();				// and re-open them at the new location
				}

				// Compute pre-seeded nucleation events.
				seedsLeft += static_cast<int> (0.499 + p.preSeededSeedDensity * geometry->area) ;
				// Note that if forbiddenZones==1, a fraction of these events will be 'non-events'


				// perform further initialization
				mtLengthHistogram.setBins(p.hiresLengthHistogramBins);
				mtLifetimeHistogram.setBins(p.hiresLifetimeHistogramBins);
				segAngleLengthHistogram.setBins(p.loresAngleHistogramBins, p.loresLengthHistogramBins);
				segAngleLifetimeHistogram.setBins(p.loresAngleHistogramBins, p.loresLifetimeHistogramBins);

				p.writeToFile();

				// if snapshot recording is enabled...
				if (p.movieEnabled)
				{
					// ... and no snapshot event has been scheduled
					if (nextSnapshotEventTime < ZERO_CUTOFF)
					{
						timeQueue.pushGlobal(p.movieFrameInterval, snapshot);
						nextSnapshotEventTime = p.movieFrameInterval;
					}
				}
				else
					nextSnapshotEventTime = 0;

				if (p.newParameterReadInterval > ZERO_CUTOFF)
				{
					timeQueue.pushGlobal(p.newParameterReadInterval, parameter_change);
					nextParameterEventTime = p.newParameterReadInterval;
				}

			}
			else
			{
				cerr << "Failed to read parameters from file. Exiting.\n";
				emergencyBreak();
				exit(-1);
			}

			timeQueue.pushGlobal(p.newParameterReadInterval, parameter_change);
			nextParameterEventTime = p.newParameterReadInterval;
			break;
	}

	return;
}




EventDescriptorIndex System::registerEventDescriptor(EventDescriptor* ei)
{
	// WARNING: maybe use a pool of numbers instead of an increasing number....this could wrap around
	EventDescriptorIndex index = EventDescriptorID++;
	EventDescriptorMap[index] = ei;

	if (EventDescriptorID == numeric_limits<EventDescriptorIndex>::max() - 10)
		timeQueue.pushGlobal(systemTime, stop);

	return index;
}

void System::unregisterEventDescriptor(EventDescriptorIndex ei)
{
	EventDescriptorMap.erase(EventDescriptorMap.find(ei));
	return;
}

EventDescriptor* System::getEventDescriptor(EventDescriptorIndex idx)
{
	map<EventDescriptorIndex, EventDescriptor*>::iterator itr = EventDescriptorMap.find(idx);
	if (itr == EventDescriptorMap.end())
		return NULL;
	else
		return itr->second;
}

void System::removeOccupiedIntersection(Intersection& is)
{
	//	cout << "removing OccupiedIntersection; Size= " << OccupiedIntersectionList.size() << "; index= " << is.occupiedListPtr->index <<"\n";
	OccupiedIntersectionList.RemoveElement(is.occupiedListPtr->index);
	is.occupiedListPtr = NULL; 
	is.mirror->second.occupiedListPtr = NULL;
	return;
}

void System::addOccupiedIntersection(IntersectionItr is)
{
	//	cout << "adding OccupiedIntersection\n";
	is->second.occupiedListPtr = OccupiedIntersectionList.create(is);
	is->second.mirror->second.occupiedListPtr = is->second.occupiedListPtr;
	return;
} 

void System::makeBinomialTable(void)
{
	int i,j;
	double temp;
	double factorial[MAXBINOM];
	temp=1;
	factorial[0] = 1;
	for (i=1 ; i<MAXBINOM ; i++)
	{
		temp *= i;
		factorial[i] = temp;
	}
	for (i=0 ; i<MAXBINOM ; i++)
		for (j=0 ; j<=i ; j++)
			binomialTable[i][j] = factorial[i]/(factorial[j]*factorial[i-j]);
	return;
}
