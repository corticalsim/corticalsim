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
    totalNucleationCount(0),
		totalUnboundNucleationCount(0),
		totalMTbasedNucleationCount(0),

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
      geometry = new Grid(Coord2D(p.geomParam1, p.geomParam2), static_cast<int>(p.geomParam3+0.4999), p.gridAspectratio, this);
      break;
    case g_wormhole:
      geometry = new Wormhole(Coord2D(p.geomParam1, p.geomParam2), this);
      break;
    case g_cylinder:
      geometry = new Cylinder(p.geomParam1, p.geomParam2, this);
      break;
    case g_gridcylinder:
      geometry = new GridCylinder(p.geomParam1, p.geomParam2, static_cast<int>(p.geomParam3+0.4999), p.gridAspectratio, this);
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

  // initialize derived parameters for double saturating nucleation
  if (p.useDoubleNucleationSaturation)
  {
    initDoubleSatPars();
  }

  // perform further initialization
  mtLengthHistogram.setBins(p.hiresLengthHistogramBins);
  mtLifetimeHistogram.setBins(p.hiresLifetimeHistogramBins);
  segAngleLengthHistogram.setBins(p.loresAngleHistogramBins, p.loresLengthHistogramBins);
  segAngleLifetimeHistogram.setBins(p.loresAngleHistogramBins, p.loresLifetimeHistogramBins);
  if (p.spatialHistogramType == sh_x || p.spatialHistogramType == sh_xy)
  {
      histX = new OneDSpatialMeasurement( p.spatialHistogramBinsX,p.spatialHistogramCountCaps,o_x,p.geomParam1,geometry->getAreaForSpatialHist(p.geomParam1,p.geomParam2),p.spatialHistogramWriteRaw);
  }
  if (p.spatialHistogramType == sh_y || p.spatialHistogramType == sh_xy)
  {
      histY = new OneDSpatialMeasurement( p.spatialHistogramBinsY,p.spatialHistogramCountCaps,o_y,p.geomParam2,geometry->getAreaForSpatialHist(p.geomParam1,p.geomParam2),p.spatialHistogramWriteRaw);
  }


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
  if (p.spatialHistogramType == sh_x || p.spatialHistogramType == sh_xy)
  {
      delete histX ;
  }
  if (p.spatialHistogramType == sh_y || p.spatialHistogramType == sh_xy)
  {
      delete histY ;
  }


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
		if (abs(geometry->regions[ridx]->totalLength - regionLength) > ZERO_CUTOFF*max(100., regionLength))
		{
			cerr << "Unacceptable drift in region length (drift = " << abs(geometry->regions[ridx]->totalLength - regionLength) << "). Exiting. [try lowering flush interval]\n";
#ifndef SLEEZY
			emergencyBreak();
			exit(-1);
#endif
		}
		geometry->regions[ridx]->totalLength = regionLength;
		systemLength += regionLength;
	}
	if (abs(totalLength - systemLength) > ZERO_CUTOFF*max(100., totalLength))//max(1.,totalLength))
	{
		cerr << "Unacceptable drift in system length (drift = " << abs(totalLength - systemLength) << "). Exiting. [try lowering flush interval]\n";
#ifndef SLEEZY
		emergencyBreak();
		exit(-1);
#endif
	}
	totalLength = systemLength;

	// reset master timer
	systemTimeOffset += timeOffset;
	systemTime = 0;

#ifdef DBG_VELOCITY
//show entire vPlusQueue
DeterministicEvent te; 
cout << "QQQ: showing vPlusQueue (flushing)\n";
while ( ! vPlusQueue.queue.empty()) {
te = vPlusQueue.pop();
  cout << "QQQ: "<< te.eventTimeDist << "\t" <<  te.infoIdx << "\t" <<  te.tag << "\t" <<  te.global_type << "\t" ; 
		if (te.infoIdx != -1)
		{
			map<EventDescriptorIndex, EventDescriptor*>::iterator eventItr = EventDescriptorMap.find(te.infoIdx);
			// Invalidated events are not removed from the queue, but the nextEvent for the segment *is* updated.
			// Therefore, check whether the tip still exists and whether the tag still matches
			if ((eventItr != EventDescriptorMap.end()) 
				&& (eventItr->second->tag == te.tag)) 
			{
				cout <<  eventItr->second->type << "\n";
			}
			else // no longer a valid event
			{
				cout << "INVALID\n";
			}
		}
    else 
    cout << "global\n";
}

cout << "QQT: showing timeQueue (flushing)\n";
while ( ! timeQueue.queue.empty()) {
te = timeQueue.pop();
  cout << "QQT: "<< te.eventTimeDist << "\t" <<  te.infoIdx << "\t" <<  te.tag << "\t" <<  te.global_type << "\t" ;
		if (te.infoIdx != -1)
		{
			map<EventDescriptorIndex, EventDescriptor*>::iterator eventItr = EventDescriptorMap.find(te.infoIdx);
			// Invalidated events are not removed from the queue, but the nextEvent for the segment *is* updated.
			// Therefore, check whether the tip still exists and whether the tag still matches
			if ((eventItr != EventDescriptorMap.end()) 
				&& (eventItr->second->tag == te.tag)) 
			{
				cout <<  eventItr->second->type << "\n";
			}
			else // no longer a valid event
			{
				cout << "INVALID\n";
			}
		}
    else 
    cout << "global\n";
}
#endif

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
	// if the pool is infinite, act as the identity function
	if (p.restrictedPool == 0)
		return dist;

	double Lmax = (p.poolDensity*geometry->area);
	double invLMax = 1.0/(p.poolDensity*geometry->area);
	double alpha = (p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size();

	if  (growing_mts.size() == 0)
	{
		if (shrinking_mts.size() == 0)
			return dist/(1.0 - totalLength*invLMax);
		else
			return ((Lmax - totalLength) - sqrt(pow(Lmax - totalLength,2) - 2.*alpha*Lmax*dist))/alpha;
	}

	double beta = p.vPlus*growing_mts.size()*invLMax;

	if ((shrinking_mts.size() == 0) && ((p.vTM/p.vPlus) < ZERO_CUTOFF))
		return -log(1.0 - growing_mts.size()*p.vPlus*dist/(Lmax-totalLength))/beta;

	double pv = (alpha/beta - totalLength)/(Lmax - alpha/beta);
	double qv = growing_mts.size()*p.vPlus*dist/(Lmax - alpha/beta);

	// a heuristic criterion for calling the LambertW function
	if (pv - qv < 250)
		return (qv-pv + LambertW(pv*exp(pv-qv)))/beta;

	// if not, do asymptotic approximation

	double L1 = pv-qv + log(pv);
	double L2 = log(L1);
	return (1.0/beta)*(qv - pv + L1 - L2 + L2/L1 + L2*(L2-2)/(2*pow(L1,2))+ L2*(6-9*L2+2*pow(L2,2))/(6*pow(L1,3))\
		+ L2*(-12+36*L2-22*pow(L2,2)+3*pow(L2,3))/(12*pow(L1,4)));

}

double System::timeToVPlus(double time)
{
	// if the pool is infinite, act as the identity function
  //if (randomGen() < 0.0001) 
    //cout << "blah " << p.vPlus << "\n";
	if (p.restrictedPool == 0)
		return time;

	double invLMax = 1.0/(p.poolDensity*geometry->area);
	double alpha = (p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size();

	if  (growing_mts.size() == 0)
		return ((1.0 - totalLength*invLMax) - 0.5*alpha*invLMax*time)*time;

	double beta = p.vPlus*growing_mts.size()*invLMax;
	double lambda_inf = invLMax*alpha/beta;
	double lambda_0 = totalLength*invLMax;

	return (1.0 - lambda_inf)*time - (lambda_0 - lambda_inf)*(1.0 - exp(-beta*time))/beta;
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
		case extraRescue:
			handleRescueEvent(true);
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
#if defined (VAR_CAT) || defined (BAND_CAT)
    double pLoc, yPos,xPos;
#endif

		// determine whether we're looking for a 'normal' or 'special' tip
//		if (randomGen.rand() < (double)growingTipsNormal/((double)growingTipsNormal + (double)p.catastropheMultiplier*growingTipsSpecial)) // added extra test for boundary case rand == 1
    if (randomGen.rand() < (double)growingTipsNormal/((double)growingTipsNormal + (double)p.catastropheMultiplier*growingTipsSpecial) or growingTipsSpecial == 0)
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


#ifdef VAR_CAT
    // Use that the system is markovian. Spontaneous catastrophe events are generated with kCatMax and discarded if rand() >= kCatNow/kCatMax. This procedure conserves both local and average rates. 
    pLoc = 1.;
    if (paramType == r_param_normal ) {
      switch (p.kCatOrient) {
        case o_x: // finding xPos and yPos could become member functions. 
          xPos = (*tipTag)->trajectory->base.x + cos((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
          xPos = geometry->xPosGridToDomain(xPos, ridx);
          pLoc = (1. + p.kCatAlpha *(sin(p.invKCatPeriod *xPos) ))/ (1. + p.kCatAlpha);
          //cout << xPos << " " << pLoc << '\n' ; 
          break;
        case o_y:
          yPos = (*tipTag)->trajectory->base.y + sin((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
          yPos = geometry->yPosGridToDomain(yPos, ridx);
          pLoc = (1. + p.kCatAlpha *(sin(p.invKCatPeriod *yPos) ))/ (1. + p.kCatAlpha);
          break;
        case o_z:
          cerr << "VAR_CAT Not implemented.\n" ;
          exit (-4);
          break;
          pLoc = (1. + p.kCatAlpha *(sin(p.invKCatPeriod *yPos) ))/ (1. + p.kCatAlpha);
          //cout << yPos << " " << pLoc << '\n' ; 
      }

    } 
    if (randomGen.rand() >= pLoc){
      return;
    } 
#endif

#ifdef BAND_CAT
    // Duplicated the above block on purpose to prevent nested ifdefs etc. 
    // Use that the system is markovian. Spontaneous catastrophe events are generated with kCatMax and discarded if rand() >= kCatNow/kCatMax. This procedure conserves both local and average rates. 
    // Caps are treated as gaps (r_param_modified)
    if (paramType == r_param_normal ) {
      double posInBand;
        if ((p.catPattern == p_band && p.kCatOrient==o_x && p.nSpirals==0) || p.catPattern == p_singleBand ) {yPos = 0;}
        else {
          yPos = (*tipTag)->trajectory->base.y + sin((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
          yPos = geometry->yPosGridToDomain(yPos, ridx);
        }
        if (p.catPattern == p_band && p.kCatOrient==o_y && p.nSpirals==0 ) {xPos = 0;}
        else {
          xPos = (*tipTag)->trajectory->base.x + cos((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
          xPos = geometry->xPosGridToDomain(xPos, ridx);
        }
        //cout << "TESTING" << "\t" << xPos << "\t" << yPos << "\n"; 
      switch (p.catPattern) {
        case p_cross:
          //if (abs (xPos - p.bla*yPos) < p.blabla2 && abs(yPos - p.bla2*xPos) < p.blabla)
                //pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
        case p_bar:
          //if (abs (xPos - p.bla*yPos) < p.blabla && abs(yPos - p.bla2*xPos) < p.blabla2)
                //pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
          break;  
        case p_singleBand:
          switch (p.kCatOrient) {
            case o_x:
              if (abs(xPos) < 0.5*p.bandGapWidth[0]) {
                pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
                //cout << "x band:  " << xPos << " pLoc: " << pLoc <<"\n";
              }
              else{
                pLoc = p.kCat[1]/p.catMax;  // gap; explictly don't assume that p.kCat[gap] > p.kCat[band]
                //cout << "x gap:  " << xPos << " pLoc: " << pLoc <<"\n";
              }
              break;
          default:
            cerr << "singleBand only implemented for orientation x\n";
            break;
         } 
        break;
      case p_band:
      switch (p.kCatOrient) {
        case o_x:
          if ( p.nSpirals > 0 ) {
            //switch (p.geometry) {
              //case (g_periodic):
                //wrapLength =  static_cast<Periodic*>(geometry)->size.x;
                //break;
              //case (g_cylinder):
                //wrapLength =  static_cast<Cylinder*>(geometry)->radius*2.*PI;
                //break;
            //}
            yPos = (*tipTag)->trajectory->base.y + sin((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
            yPos = geometry->yPosGridToDomain(yPos, ridx);
          }
          else{
            yPos=0;
          }
          xPos = (*tipTag)->trajectory->base.x + cos((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
          xPos = geometry->xPosGridToDomain(xPos, ridx);
          //if ( abs(0.5*(p.bandGapWidth[0]+p.bandGapWidth[1]) - fmod(abs(xPos+p.spiralPitch*yPos),(p.bandGapWidth[0]+p.bandGapWidth[1]))) <= 0.5*p.bandGapWidth[0]) {
            //pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
          //}
          //else{
            //pLoc = p.kCat[1]/p.catMax;  // gap; explictly don't assume that p.kCat[gap] > p.kCat[band]
          //}
          pLoc = p.kCat[1]/p.catMax;  // gap; explictly don't assume that p.kCat[gap] > p.kCat[band]
          if ( p.nSpirals != 0 ) {
            for (int k= 0 ; k < p.nSpirals ; k++ ) {
              //if ( abs(0.5*(p.bandGapWidth[0]+p.bandGapWidth[1]) - fmod(abs(xPos+k*wrapLength/p.nSpirals+p.spiralPitch*yPos),(p.bandGapWidth[0]+p.bandGapWidth[1]))) <= 0.5*p.bandGapWidth[0]) (
              posInBand = fmod(abs(-(0.5 + k)*p.projectedPeriod + xPos - yPos*p.nSpirals*p.projectedPeriod/p.wrapLength),p.projectedPeriod) ;
              if ( posInBand <= 0.5*p.projectedBand || posInBand >= p.projectedPeriod - 0.5*p.projectedBand ) {
                pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
                //cout << "SPIRAL " << p.nSpirals << "\tx: " << xPos << "\ty: " << yPos << "\n";
              }
            }
          }
          else {
            posInBand = fmod(abs(-0.5*p.projectedPeriod + xPos - yPos*p.nSpirals*p.projectedPeriod/p.wrapLength),p.projectedPeriod) ; 
            if ( posInBand <= 0.5 * p.projectedBand || posInBand >= p.projectedPeriod - 0.5*p.projectedBand ) {
              //if ( abs(0.5*(p.bandGapWidth[0]+p.bandGapWidth[1]) - fmod(abs(xPos+p.spiralPitch*yPos),(p.bandGapWidth[0]+p.bandGapWidth[1]))) <= 0.5*p.bandGapWidth[0]) (
              pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
            }
          }
          //cout << xPos << " " << pLoc << '\n' ; 
          break;
        case o_y:
          if ( p.nSpirals > 0 ) {
            xPos = (*tipTag)->trajectory->base.x + cos((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
            xPos = geometry->xPosGridToDomain(xPos, ridx);
          }
          else{
            xPos=0;
          }
          yPos = (*tipTag)->trajectory->base.y + sin((*tipTag)->trajectory->base.angle) *(*tipTag)->position();
          yPos = geometry->yPosGridToDomain(yPos, ridx);
          pLoc = p.kCat[1]/p.catMax;  // gap; explictly don't assume that p.kCat[gap] > p.kCat[band]
          if ( p.nSpirals != 0 ) {
            for (int k= 0 ; k < p.nSpirals ; k++ ) {
              posInBand = fmod(abs(-(0.5 + k)*p.projectedPeriod + yPos - xPos*p.nSpirals*p.projectedPeriod/p.wrapLength),p.projectedPeriod); 
              if ( posInBand <= 0.5*p.projectedBand || posInBand >= p.projectedPeriod - 0.5*p.projectedBand ) {
                pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
                //cout << "SPIRAL " << p.nSpirals << "\tx: " << xPos << "\ty: " << yPos << "\n";
              }
            }
          }
          else {
            posInBand= fmod(abs(-0.5*p.projectedPeriod + yPos ),p.projectedPeriod);
            if ( posInBand <= 0.5*p.projectedBand || posInBand >= p.projectedPeriod - 0.5*p.projectedBand ) {
              pLoc = p.kCat[0]/p.catMax;  // band; explictly don't assume that p.kCat[gap] > p.kCat[band]
            }
          }
          //if ( abs(0.5*(p.bandGapWidth[0]+p.bandGapWidth[1]) - fmod(abs(yPos+p.spiralPitch*xPos),(p.bandGapWidth[0]+p.bandGapWidth[1]))) <= 0.5*p.bandGapWidth[0]) {
            //pLoc = p.kCat[0]/p.catMax; //band
          //}
          //else{
            //pLoc = p.kCat[1]/p.catMax; // gap
          //}
          break;
        case o_z:
          cerr << "BAND_CAT Not implemented.\n" ;
          exit (-4);
          break;
          //cout << yPos << " " << pLoc << '\n' ; 
      }
      break;
      }

    } 
    else {
      pLoc = p.kCat[1]/p.catMax; // Treat caps as "gap"
    }
    if (randomGen.rand() >= pLoc){
      return;
    } 
#endif

    (*tipTag)->mt->catastrophe();

    return;
}



#ifndef NO_INLINE
inline 
#endif
void System::handleRescueEvent(bool extra)
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

		// ignore fraction of the extra events depending on function
		if (extra) {
			double theta = (*tipTag) -> trajectory -> base.angle; // this implementation cannot account for polarity. 
			switch (p.extraRescueFunction) {
				case r_cos: // 1/2 + 1/2 *cos(2 (theta - kResExtraAngle))
					if (randomGen.rand() > 0.5 + 0.5*cos(2*(theta - p.kResExtraAngle)))	
					{return;}
					break;
				default:
					break;
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

void System::randomPositionOnMicrotubuleInRegion(int ridx, double cutLength, double& randomPos, Segment*& randomSeg)
{
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


void System::randomSegmentAtMetaIntersection(double posOnTrajectory, Trajectory* tr, int occupancy, Segment*& randomSeg)///////////////////////////////////////////////////////////////////////////////////////////////////
{

	vector<bool> segAtIntersection;
	int count(0);
	list<Segment*>::iterator seg = tr->segments.begin();
	while ( seg != tr->segments.end() )                                                              // here we tag all segments on the trajectory that cross the MetaIntersection
	{

		(**seg).mt->updateLength();
		if ((**seg).crossesPoint(posOnTrajectory))
		{
			segAtIntersection.push_back(true);
			++count;
		}
		else
		{
			segAtIntersection.push_back(false);
		}
		++seg;
	}


	if ( count != occupancy )
	{
		randomSeg = *(tr->segments.begin());
		
		cerr << "Error in selecting a random segment that crosses a meta intersection. Detected: " << count << " against an occupancy number: " << occupancy << ". Abort! "<< tr->segments.size() << "\n";
		exit(2);
	}
	
	int RVseg = randomGen.randDblExc(count);//randomGen.randDblExc(occupancy);                                                       // here we select one random segment among those that have been tagged
	seg = tr->segments.begin();
	count=0;    // <---------------------------------------------------------------------------------------------------- check this out, probably correct
	for (int i=0; i<segAtIntersection.size(); ++i)
	{
		if ( segAtIntersection.at(i) )
		{
			if ( RVseg==count )
				break;
			else
				++count;
		}
		++seg;
	}
	
	randomSeg = *seg;
	return;
}


void System::selectRandomRegionProportionalBy(string quantitytype, double quantitymax, int& ridx, double& cutLength)
{
  ridx = 0;
  double regionalquantity;

	// select position at which severing will take place
	cutLength = randomGen.randDblExc(quantitymax);


	// first, select relevant region. Go from both sides to reduce time.
	if (cutLength < 0.5*quantitymax)
	{
		ridx = 0;
		while (true)
		{
      if (quantitytype == "length")
      {
  			geometry->regions[ridx]->updateRegionLength();
        regionalquantity = geometry->regions[ridx]->totalLength;
      }
      else if (quantitytype == "area")
      {
        regionalquantity = geometry->regions[ridx]->area;
      }
      else {exit(-1);}
			if (cutLength < regionalquantity)
				break;
			cutLength -= regionalquantity;
			ridx++;
			if (ridx == geometry->regions.size())
			{
				// apparently, cutLength was *just* too long...
				cutLength -= ZERO_CUTOFF;
#ifdef DBG_ACID_TEST
				if (cutLength > regionalquantity)
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
		cutLength = quantitymax - cutLength;
		while (true)
		{
      if (quantitytype == "length")
      {
  			geometry->regions[ridx]->updateRegionLength();
        regionalquantity = geometry->regions[ridx]->totalLength;
      }
      else if (quantitytype == "area")
      {
        regionalquantity = geometry->regions[ridx]->area;
      }
      else {exit(-1);}
			if (cutLength < regionalquantity)
				break;
			cutLength -= regionalquantity;
			ridx--;
			if (ridx == -1)
			{
				// apparently, cutLength was *just* too long...
				cutLength -= ZERO_CUTOFF;
#ifdef DBG_ACID_TEST
				if (cutLength > regionalquantity)
				{
					cerr << "ERROR: cannot locate intersection location\n";
					exit(-1);
				}
#endif
				break;
			}
		}
		//and wrap it back
		cutLength = regionalquantity - cutLength;
	}
}

void System::randomPositionOnMicrotubule(double& randomPos, Segment*& randomSeg)
{
	int ridx;
  double cutLength;

  // Find a random region proportional by microtubule length within the region (defines ridx and cutLength)
	selectRandomRegionProportionalBy("length", totalLength, ridx, cutLength);

	// Now that the area has been selected, cycle over trajectories and their segments
  randomPositionOnMicrotubuleInRegion(ridx, cutLength, randomPos, randomSeg);
	return;
}


bool System::handleRegionalNucleationSaturation(int ridx, double& randomPos, Segment*& randomSeg)
{
  // Handle nucleation in a random region proportional to region area
  // Uses same region index (ridx) assigned for case of isotropic nucleation

  // Update region length and determine local density
  geometry->regions[ridx]->updateRegionLength();
  double regionaldensity = geometry->regions[ridx]->totalLength / geometry->regions[ridx]->area;

  // Calculate regional bound nucleation rate as a fraction of the max bound nucleation rate
  double frn = regionaldensity / (p.nucleationHalfIsotropicDensity + regionaldensity);

  // Discard overbooked nucleations
  if (randomGen.randExc() >= frn)
    return false;

  // Determine segment and position for bound nucleation
  double cutLength = randomGen.randDblExc(geometry->regions[ridx]->totalLength);
  randomPositionOnMicrotubuleInRegion(ridx, cutLength, randomPos, randomSeg);
  return true;
}

void System::initDoubleSatPars()
{
  double area = geometry->area;
  if (p.forbiddenZones and p.geometry == g_gridcylinder)
    area -= 2*PI*pow(static_cast<GridCylinder*>(geometry)->radius,2);
  double rnmax_target = p.kNuc/p.rn_targetratio;
  double nNCocc_target = rnmax_target*p.occupancytimeNC*area;
  p.nNCmax = static_cast<int>(p.kNuc*nNCocc_target / (p.kNuc - rnmax_target));
  double rnbase_target = p.f_unbound * rnmax_target;
  double f_occ_target = (p.nNCmax - nNCocc_target) / p.nNCmax;
  p.rn_base0 = rnbase_target/f_occ_target;
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
		if (p.crossSeveringTop || randomGen.rand() < p.crossSeveringTopFraction )//(randomGen.randInt(1)==0))
			// if MT (bundle) on top should be cut, select that one, or (if random) pick it with p.crossSeveringTopFraction chance.
			// Note that no effort is made to ensure proportional selection for bundles of different thickness.
			cutIS = cutOccIS->intersectionToCut; 
		else
			cutIS = cutOccIS->intersectionToCut->ISREF.mirror;
		double diffAngle = cutIS->ISREF.otherTrajectory->base.region->intersectionAngle(
				cutIS->ISREF.otherTrajectory,
				cutIS->ISREF.mirror->ISREF.otherTrajectory);
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
		cutSegNumber = randomGen.randInt(cutIS->ISREF.occupancy - 1)  ; // nummers 0 - occupancy -1 .
		if (cutSegNumber < cutIS->ISREF.occupancy/2)	// search from the beginning.
		// In principle, all searches could start from the beginning. Also starting from the end should save (on average) half of the time 
		{
			testSeg = cutIS->ISREF.mirror->ISREF.otherTrajectory->segments.begin();
			finalTestSeg = cutIS->ISREF.mirror->ISREF.otherTrajectory->segments.end();
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
			testSeg = cutIS->ISREF.mirror->ISREF.otherTrajectory->segments.end();
			finalTestSeg = cutIS->ISREF.mirror->ISREF.otherTrajectory->segments.begin();
			while (testSeg != finalTestSeg)
			{
				testSeg--;
				if ((**testSeg).crossesIntersection(cutIS))
				{				
					cutSegNumber++;
					if (cutSegNumber == cutIS->ISREF.occupancy)
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


bool System::handleReducedGapNucleation(SurfaceVector &sv, bool useIsotropic, Segment* nucSeg, double* pospointer)
{
// reject fraction of nucleation events if in gap region. Note that this reduces the total nucleation rate. 
  double posInBand,xPos,yPos;
  double pos;
  if (not useIsotropic)
    pos = *pospointer;

  #ifdef BAND_CAT
  if (p.nSpirals != 0) {
    cerr << "reducedGapNucleation not implemented for spirals!\n";
    exit(-16);
  }
  #endif

  int ridx;
  if (useIsotropic)
    ridx = sv.region->geometryRegionIndex ;
  else
    ridx = nucSeg->trajectory->base.region->geometryRegionIndex;
  
  #ifdef BAND_CAT
  switch (p.kCatOrient) {
    case o_x:
      if (useIsotropic)
        xPos = sv.x;
      else
        xPos = nucSeg->trajectory->base.x + cos(nucSeg->trajectory->base.angle) *pos;
      xPos = geometry->xPosGridToDomain(xPos, ridx);

      switch (p.catPattern) {
        case p_band:
          posInBand = fmod(abs(-0.5*p.projectedPeriod + xPos ),p.projectedPeriod) ; 
          if ( posInBand > 0.5 * p.projectedBand && posInBand < p.projectedPeriod - 0.5*p.projectedBand && randomGen.randDblExc() > p.gapNucleationAcceptFraction ) {
            return true; // reject gap with certain probability.
          }
          break;
        case p_singleBand:
          if (abs(xPos) > 0.5*p.bandGapWidth[0] && randomGen.randDblExc() > p.gapNucleationAcceptFraction ) {
            return true; // reject gap with certain probability.
          }
          break;
      }
      break;
    case o_y:
      if (useIsotropic)
        yPos = sv.y;
      else
        yPos = nucSeg->trajectory->base.y + sin(nucSeg->trajectory->base.angle) *pos;
      yPos = geometry->yPosGridToDomain(yPos, ridx);
      switch (p.catPattern) {
        case p_band:
          posInBand= fmod(abs(-0.5*p.projectedPeriod + yPos ),p.projectedPeriod);
          if ( posInBand > 0.5*p.projectedBand && posInBand < p.projectedPeriod - 0.5*p.projectedBand && randomGen.randDblExc() > p.gapNucleationAcceptFraction ) {
            return true;  // reject gap with certain probability.
          }
          break;
       case p_singleBand:
          if (abs(yPos) > 0.5*p.bandGapWidth[0] && randomGen.randDblExc() > p.gapNucleationAcceptFraction ) {
            return true; // reject gap with certain probability.
          }
          break;
      }
      break;
    default:
      cerr << "Reduced nucleation in gaps not implemented for this orientation\n";
      break;

  }
  #endif
  return false;
}


void System::handleBoundNucleation(SurfaceVector& sv, TrajectoryVector& tv, Segment* nucSeg, double& pos)
{
  double rand, rand2, preTheta, cosTheta, temp, overrideAngle;
  rand = randomGen.randDblExc();
  if (p.ellipseForwardAlongMT && rand >= p.ellipseLeftFraction + p.ellipseRightFraction)
  {
    // Shifting (redistributeNucleationOverBands) not implemented on purpose: would create bookkeeping nightmare. 
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
    // Pick a value of ellipseEpsilon for branched or (anti)parallel nucleation
    double localEpsilon;
    if ( rand < p.ellipseLeftFraction + p.ellipseRightFraction ) 
    {
      localEpsilon = p.ellipseEpsilon;
    } 
    else 
    {
      localEpsilon = p.ellipseEpsilonAlongMT;
    }
    sv = nucSeg->trajectory->base;
    sv.region->translateVector(sv, pos);
    preTheta = randomGen.randDblExc()*2.*PI;
    cosTheta = cos(preTheta);
    rand2 = randomGen.randDblExc();
    rand2 = sqrt(rand2);
    temp = rand2*sqrt((1.-localEpsilon*localEpsilon)*(1.-cosTheta*cosTheta));
    if (preTheta > PI)
    {
      temp = - temp;
    } 
    overrideAngle = atan2(temp, (localEpsilon + rand2*cosTheta));

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

#ifdef BAND_CAT 
    if (p.redistributeNucleationOverBands){
      // IMPORTANT: SHIFTING NOW ONLY ON OFF-(PARENT)TRAJECTORY EVENTS. Parameter check: does not run if p.redistributeNucleationOverBands && p.ellipseForwardAlongMT. 
      if (sv.region->parameterType  == r_param_normal ) { 
        int shift =  randomGen.randInt(p.redistributeShiftNumber -1);
        //cout << "Testing: shift: " << shift << " sv.x: " << sv.x << " sv.y: " << sv.y << "\n";
        if (shift > 0) {
          switch (p.kCatOrient) {
            case o_x:
              if ( p.geometry == g_grid || p.geometry == g_gridcylinder ) {
                geometry->shiftSurfaceVectorNoWrap(sv,shift*p.redistributeShift,(shift - p.redistributeShiftNumber) * p.redistributeShift, 0, 0); 
              } 
              else {
                if (sv.x + shift*p.redistributeShift > p.redistributeShiftMax)
                  sv.x += (shift - p.redistributeShiftNumber) * p.redistributeShift;
                else
                  sv.x += shift*p.redistributeShift;
              }
              break; 
            case o_y:
              if ( p.geometry == g_grid || p.geometry == g_gridcylinder ) {
                geometry->shiftSurfaceVectorNoWrap(sv,0,0,shift*p.redistributeShift,(shift - p.redistributeShiftNumber) * p.redistributeShift); 
              } 
              else {
                if (sv.y + shift*p.redistributeShift > p.redistributeShiftMax)
                  sv.y += (shift - p.redistributeShiftNumber) * p.redistributeShift;
                else
                  sv.y += shift*p.redistributeShift;
              } 
              break;
          }
        }
        //cout << "Testing again: shift: " << shift << " sv.x: " << sv.x << " sv.y: " << sv.y << " region: " << sv.region->geometryRegionIndex << "\n";
        //cout << "Testing again: shift: " << shift << " sv.x: " << sv.x << " sv.y: " << sv.y << " region: " << sv.region->geometryRegionIndex << "\n";
      }
    }
#endif

    sv.region->translateVector(sv, ZERO_CUTOFF);
    tv = geometry->createTrajectory(sv);
  }
  if ( true ) // extreem smerige constructie om gezeik van de compiler te voorkomen :-((
  {
    Microtubule* mt = growing_mts.create(this, tv);
		
    if (p.outputNucPos == "X" or p.outputNucPos == "XY")
    {
      double xPosMeas;
      xPosMeas = sv.x;
      xPosMeas = geometry->xPosGridToDomain(xPosMeas, sv.region->geometryRegionIndex);
      nucleationXpositions.push_back(xPosMeas);
    }
    if (p.outputNucPos == "Y" or p.outputNucPos == "XY")
    {
      double yPosMeas;
      yPosMeas = sv.y;
      yPosMeas = geometry->yPosGridToDomain(yPosMeas, sv.region->geometryRegionIndex);
      nucleationYpositions.push_back(yPosMeas);
    }
  }
  return;
}


double System::handleBiasedNucleationAngle(SurfaceVector& sv, double overrideAngle, bool preSeeded)
{
  if (p.nucleationBiasType == nbias_cross ) {
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
      return overrideAngle;
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
  }
  else {  // biasType: ellipsePolar or ellipseApolar
    // Distribution similar to ellipse distributions used in nuc_ellipse, but uses parameter "nucleationAlpha" in stead of "ellipseEpsilon" for full control. 
    int rot=0;
    double preTheta = randomGen.randDblExc()*2.*PI;
    double cosTheta = cos(preTheta);
    double rand2 = randomGen.randExc();  // use randExc in stead of randDblExc
    if (p.nucleationBiasType == nbias_ellipseApolar ) {  // This saves generating 1 random number; does not introduce bias (only reduces resolution of random angles by factor 2.)
      rand2 *= 2;
      if (rand2 >= 1) {
        rand2 -= 1; 
        rot=1;
      }
    } 
    rand2 = sqrt(rand2);
    double temp = rand2*sqrt((1.-p.nucleationAlpha*p.nucleationAlpha)*(1.-cosTheta*cosTheta)); 
    //temp = rand2*sqrt((1.-p.ellipseEpsilon*p.ellipseEpsilon)*(1.-cosTheta*cosTheta));
    if (preTheta > PI)
    {
      temp = - temp;
    } 
    overrideAngle = atan2(temp, (p.nucleationAlpha + rand2*cosTheta)) + rot*PI + p.nucleationBiasAngle;
    //overrideAngle = atan2(temp, (p.ellipseEpsilon + rand2*cosTheta)) + rot*PI;
    // for safety: bring back to interval [0..2 PI]
    while (overrideAngle > 2*PI)
      overrideAngle -= 2*PI;
    while (overrideAngle < 0)
      overrideAngle += 2*PI;
  }
  return overrideAngle;
}


#ifndef NO_INLINE
inline 
#endif
void System::handleNucleationEvent(bool preSeeded)
{
/* exit immediately if there is a finite tubulin pool and it is exhausted to the point where vPlus is essentially the same as vTM. 
	Nucleating in this case would lead to moving microtubules with zero length, causing ordering problems on collisions 
	[the trailing end can sometimes collide first].

  */

  if ((p.restrictedPool != 0) && (p.vTM/(p.vPlus*(1. - totalLength/(p.poolDensity*geometry->area))) > 0.9999))
    return;

  TrajectoryVector tv;
  SurfaceVector sv = geometry->randomSurfaceVector();
  Segment* nucSeg;
  NucleationType thisNuc;
  double pos,posOnSeg;
#ifdef BAND_CAT
  double posInBand,xPos,yPos;
#endif
  double overrideAngle;
  bool boundnucleation = false;

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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// FROM HERE
  switch (thisNuc)
  {
    case nuc_isotropic:
      overrideAngle = sv.angle;
#ifdef BAND_CAT
      if (p.reducedGapNucleation) {
        if (handleReducedGapNucleation(sv, true, nullptr, nullptr))
          {
            return; // reject gap nucleation with certain probability
          }
      }
#endif
totalUnboundNucleationCount++;
      break;
    case nuc_discreteAngles:

      overrideAngle = p.nucleationAngles[randomGen.randInt(p.discreteAngleNumber-1)] + PI*randomGen.randInt(1);

      break;
      
    case nuc_aster:
    {  
	  if ( p.rvRejectionUnbound == 1 && systemTime+systemTimeOffset < 50 )           // This is to make sure that, if only Aster nucleations are to occur, we have few initial MTs where the MT-based nucleation can start off
	  {
	  	//cout << systemTime+systemTimeOffset << endl;
	  	overrideAngle = sv.angle;
		#ifdef BAND_CAT
		if (p.reducedGapNucleation) 
		{
		  if (handleReducedGapNucleation(sv, true, nullptr, nullptr))
		  {
		      return; // reject gap nucleation with certain probability
		  }
		}
		#endif
		totalUnboundNucleationCount++;
		if (p.useDoubleNucleationSaturation and (not preSeeded))
		  occupiedNCs.push_back(systemTime + systemTimeOffset);
		break;
	  }

	  if (p.useDoubleNucleationSaturation)
	  {
	  	while ((not occupiedNCs.empty()) and occupiedNCs.front() < systemTime + systemTimeOffset - p.occupancytimeNC)
		  occupiedNCs.pop_front();
	  	int nNCocc = occupiedNCs.size();
	  	double current_fNCfree = (p.nNCmax - nNCocc) / static_cast<double>(p.nNCmax);
	  	if ( randomGen.randDblExc() >= current_fNCfree )
		  return; // Discard overbooked nucleations
	  }
	  
  	  Aster aster = geometry->createAster(p.numberOfRays, sv);                                 // just selects the angles;  

	  vector<double> findLatticeProbabilities;
	  double findLatticeNorm(0.);                                          // normalization constant for the probability for nuc complex to attach to any MT lattice
	
	  for (int i=0; i<p.numberOfRays; ++i)
	  {
	  	SurfaceVector sv2 = sv;        // check if is shallow or deep copy            --------------- might be: angles internally are betwee -p/2 and p/2
	  	sv2.angle = aster.angles.at(i);
	  	MetaTrajectory* mTraj = geometry->metaTrajectories.create(this, p.maxMetaTrajectoryLength, sv2);
	  	if ( mTraj->metaTrajectoryLength <= p.maxMetaTrajectoryLength )                                               // this isn't very clean
	  	{
		  findLatticeProbabilities.push_back( exp( -p.unboundNucleationRate*pow( mTraj->metaTrajectoryLength, 2 ) / (2*p.diffusionCoefficient) ) );  //this is the probability to reach the lattice before nucleation occurs
	  	  findLatticeNorm += exp( -p.unboundNucleationRate*pow( mTraj->metaTrajectoryLength, 2 ) / (2*p.diffusionCoefficient) );
	  	}
	  	else
	  	{
	  		findLatticeProbabilities.push_back(0);
	  	}
	  	mTraj->removeTrajectories();    // because we don't need them anymore
	  }
	  if ( p.numberOfRays != 0 )
	  	findLatticeProbabilities.push_back( p.numberOfRays - findLatticeNorm );
 
  
	  double findLatticeCumulative(0.);
	  int rayIdentifier(0);
	  double rvNucDir = randomGen.randDblExc();
	  for ( vector<double> ::iterator rvm = findLatticeProbabilities.begin(); rvm != findLatticeProbabilities.end(); ++rvm )
	  {
	  	findLatticeCumulative += *rvm;
	  	if ( rvNucDir <= findLatticeCumulative / p.numberOfRays )
	  		break;
	  	++rayIdentifier;
	
	  }
	  
//	  cout << rayIdentifier << "	" << endl;
	  	
	  if ( rayIdentifier == p.numberOfRays + 1 )
	  {
	  	cerr << "Possible error in the selection of nucleation type/position, rejecting.\n";
	  	geometry->metaTrajectories.removeAll();
	  	return;
	  }
	  if ( rayIdentifier == p.numberOfRays )    // THAT'S CORRECT, CHECKED 3 BILLION TIMES, LEAVE IT AS IT IS
	  {
	  	double RVnucUnbound = randomGen.randDblExc();
	  	if ( RVnucUnbound > p.rvRejectionUnbound ) 
	    	{	

			overrideAngle = sv.angle;

			#ifdef BAND_CAT
			if (p.reducedGapNucleation) 
			{
			  if (handleReducedGapNucleation(sv, true, nullptr, nullptr))
			  {
			      return; // reject gap nucleation with certain probability
			  }
			}
			#endif

			totalUnboundNucleationCount++;
			if (p.useDoubleNucleationSaturation and (not preSeeded))
			  occupiedNCs.push_back(systemTime + systemTimeOffset);
	    		
	    		#ifdef MTBASED_NUCLEATION_PROBABILITY
			localDensity locDensity;
			MetaTrajectory* mtr2 = geometry->metaTrajectories.first();
			while ( mtr2 != NULL )
			{
				locDensity.nearbyRegions.insert( locDensity.nearbyRegions.begin(), mtr2->crossedRegions.begin(), mtr2->crossedRegions.end() );
				mtr2 = mtr2->next();
			}
			sort( locDensity.nearbyRegions.begin(), locDensity.nearbyRegions.end() );
			locDensity.nearbyRegions.erase( unique( locDensity.nearbyRegions.begin(),locDensity.nearbyRegions.end() ), locDensity.nearbyRegions.end() );
			locDensity.nearbyDensity = 0;
			for( int i=0; i < locDensity.nearbyRegions.size(); ++i )                                                 // WARNING: only works when regions have identical area
				locDensity.nearbyDensity += locDensity.nearbyRegions[i]->totalLength / ( locDensity.nearbyRegions.size() * locDensity.nearbyRegions[i]->area );      
			
			regionalDensityVectorUnbound.push_back( sv.region->totalLength / sv.region->area );
			nearbyDensityVectorUnbound.push_back( locDensity.nearbyDensity );
			globalDensityVectorUnbound.push_back( totalLength / geometry->area );
			#endif
	  	}
	  	else 
	  	{
	  		geometry->metaTrajectories.removeAll();
	  		return;
	  	}
	  }
	  else
	  {
	  	MetaTrajectory* mtr = geometry->metaTrajectories.first();
	  	int rayCounter(0);
	  	#ifdef DBG_ASTER
	  	if ( geometry->metaTrajectories.size() != p.numberOfRays )
	  		cerr << "DBG_ASTER::Number of meta trajectories and number of rays do not match.\n";
	  	#endif
	  	//if ( rayIdentifier > geometry->metaTrajectories.size() ) cerr << geometry->metaTrajectories.size() << "Come ci sono finito qui?\n"; ///////////////////////keep it, numerical error, push it back to rayId -1
	  	
//	  	if ( rayIdentifier != p.numberOfRays )             // useless, we are already in that statement
//	  	{
	  	for (int i=0; i<rayIdentifier; ++i)
	  		mtr = mtr->next();
//	  	}	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  	
		if ( mtr->doesItIntersectOccupied )
		{
			int occ(0);
			Trajectory* otherTr = mtr->nsc.motherTrajectory;
			list<Segment*>::iterator seg = otherTr->segments.begin();
			for (int i=0; i<otherTr->segments.size(); ++i)
			{
				(**seg).mt->updateLength();
				if ( (**seg).crossesPoint( mtr->nsc.posOnMother ) )
				{
					++occ;
				}
				++seg;
			}
			if ( occ == 0 )
				mtr->doesItIntersectOccupied = false;
			else
				mtr->nsc.numOfSegmentsInBundle = occ;
		}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  	if ( mtr->doesItIntersectOccupied && rayIdentifier != p.numberOfRays )
	  	{
	  		double RVnucBound = randomGen.randDblExc();
	  		if ( RVnucBound > p.rvRejectionBound )
	  		{
	  			randomSegmentAtMetaIntersection(mtr->nsc.posOnMother, mtr->nsc.motherTrajectory, mtr->nsc.numOfSegmentsInBundle, nucSeg);   //select the mother microtubule in case the meta trajectory intersects with 	a bundle

	  			pos = mtr->nsc.posOnMother;                              //should be nucSeg->start + (mtr->ncs.posOnMother - nucSeg->start);
	  		
	  			#ifdef MTBASED_NUCLEATION_PROBABILITY
				localDensity locDensity;
				MetaTrajectory* mtr2 = geometry->metaTrajectories.first();
				while ( mtr2 != NULL )
				{
					locDensity.nearbyRegions.insert( locDensity.nearbyRegions.begin(), mtr2->crossedRegions.begin(), mtr2->crossedRegions.end() );
					mtr2 = mtr2->next();
				}
				sort( locDensity.nearbyRegions.begin(), locDensity.nearbyRegions.end() );
				locDensity.nearbyRegions.erase( unique( locDensity.nearbyRegions.begin(),locDensity.nearbyRegions.end() ), locDensity.nearbyRegions.end() );
				locDensity.nearbyDensity = 0;
				for( int i=0; i < locDensity.nearbyRegions.size(); ++i )
				{                //WARNING: only works when regions have an identical area
				  locDensity.nearbyDensity += locDensity.nearbyRegions[i]->totalLength / ( locDensity.nearbyRegions.size() * locDensity.nearbyRegions[i]->area ); 

				}
				regionalDensityVector.push_back( sv.region->totalLength / sv.region->area );
				nearbyDensityVector.push_back( locDensity.nearbyDensity );
				globalDensityVector.push_back( totalLength / geometry->area );
				#endif
				
				#ifdef BAND_CAT
				if (p.reducedGapNucleation)
				{
				   if (handleReducedGapNucleation(sv, false, nucSeg, &pos))
				   {
				     return; // reject gap nucleation with certain probability // NB this now only happens to bound nucleations!
				   }
				}
				#endif
				handleBoundNucleation(sv, tv, nucSeg, pos); // create the new microtubule with the right position and orientation HERE I NEED TO PUT THE sv OF THE SEGMENT
				totalNucleationCount++;
				++totalMTbasedNucleationCount;
				if (p.useDoubleNucleationSaturation and (not preSeeded))
				  occupiedNCs.push_back(systemTime + systemTimeOffset);
				geometry->metaTrajectories.removeAll();
				return; 

			}
			else
			{
				geometry->metaTrajectories.removeAll();
	  			return;     //nucleation rejected
			}
	
	  	}
	  	else
	  	{
	  		//if (rayIdentifier == p.numberOfRays) cout << "cacca " << endl;
	  		double RVnucUnbound = randomGen.randDblExc();
	  		if ( RVnucUnbound > p.rvRejectionUnbound ) 
	  	  	{	
				overrideAngle = sv.angle;
				#ifdef BAND_CAT
				if (p.reducedGapNucleation) 
				{
				  if (handleReducedGapNucleation(sv, true, nullptr, nullptr))
				  {
				      return; // reject gap nucleation with certain probability
				  }
				}
				#endif

				totalUnboundNucleationCount++;
				if (p.useDoubleNucleationSaturation and (not preSeeded))
			  		occupiedNCs.push_back(systemTime + systemTimeOffset);
			  
	  	  		#ifdef MTBASED_NUCLEATION_PROBABILITY
				localDensity locDensity;
				MetaTrajectory* mtr2 = geometry->metaTrajectories.first();
				while ( mtr2 != NULL )
				{
					locDensity.nearbyRegions.insert( locDensity.nearbyRegions.begin(), mtr2->crossedRegions.begin(), mtr2->crossedRegions.end() );
					mtr2 = mtr2->next();
				}
				sort( locDensity.nearbyRegions.begin(), locDensity.nearbyRegions.end() );
				locDensity.nearbyRegions.erase( unique( locDensity.nearbyRegions.begin(),locDensity.nearbyRegions.end() ), locDensity.nearbyRegions.end() );
				
				locDensity.nearbyDensity = 0;
				for( int i=0; i < locDensity.nearbyRegions.size(); ++i )                                                 // WARNING: only works when regions have identical area
					locDensity.nearbyDensity += locDensity.nearbyRegions[i]->totalLength / ( locDensity.nearbyRegions.size() * locDensity.nearbyRegions[i]->area );      
			
				regionalDensityVectorUnbound.push_back( sv.region->totalLength / sv.region->area );
				nearbyDensityVectorUnbound.push_back( locDensity.nearbyDensity );
				globalDensityVectorUnbound.push_back( totalLength / geometry->area );
				#endif
	  		}
	  		else
	  		{
	  			geometry->metaTrajectories.removeAll();
	  			return;     //nucleation rejected
	  		}
	  		//thisNuc = nuc_isotropic;   QUESTO NO
	  	}
	  }
	  
	  geometry->metaTrajectories.removeAll();
          break;
    }      
    case nuc_ellipse:

      double current_f_unbound;
      if (p.useDoubleNucleationSaturation)
      {
        while ((not occupiedNCs.empty()) and occupiedNCs.front() < systemTime + systemTimeOffset - p.occupancytimeNC)
          occupiedNCs.pop_front();
        int nNCocc = occupiedNCs.size();
        double current_fNCfree = (p.nNCmax - nNCocc) / static_cast<double>(p.nNCmax);
        double current_rn_max = p.kNuc * current_fNCfree;
        double frn_max_global = current_rn_max / p.kNuc;
        if (randomGen.randDblExc() >= frn_max_global)
          return; // Discard overbooked nucleations
        double current_rn_base = p.rn_base0 * current_fNCfree;
        current_f_unbound = current_rn_base / current_rn_max;
      }
      else if (p.useRegionalNucleationSaturation)
      {
        current_f_unbound = p.f_unbound;
      }
      else
      {
        current_f_unbound = p.nucleationHalfIsotropicDensity/(totalLength/geometry->area + p.nucleationHalfIsotropicDensity);
      }
      if (randomGen.randDblExc() >= current_f_unbound)
        boundnucleation = true;

      if (boundnucleation)
      {
        if (p.useRegionalNucleationSaturation or p.useDoubleNucleationSaturation)
        {
          if (not handleRegionalNucleationSaturation(sv.region->geometryRegionIndex, posOnSeg, nucSeg)) // Get segment and position if bound nucleation is to be executed
            return; // Discard overbooked nucleations
        }
        else
        {
          randomPositionOnMicrotubule(posOnSeg, nucSeg);
        }
        pos = nucSeg->start + nucSeg->dir*posOnSeg;
#ifdef BAND_CAT
        if (p.reducedGapNucleation) {
          if (handleReducedGapNucleation(sv, false, nucSeg, &pos)){
            return; // reject gap nucleation with certain probability // NB this now only happens to bound nucleations!
          }
        }
#endif
        handleBoundNucleation(sv, tv, nucSeg, pos); // create the new microtubule with the right position and orientation
        totalNucleationCount++;
	totalMTbasedNucleationCount++;
        if (p.useDoubleNucleationSaturation and (not preSeeded))
        {
          occupiedNCs.push_back(systemTime + systemTimeOffset);
        }
        return; 
      }
      else // unbound nucleation
      {
        if ( p.ellipseReducedFreeRate == 0 || randomGen.randDblExc() < p.ellipseReducedFreeRateAcceptFraction ) {
          if ( abs(p.nucleationAlpha) < ZERO_CUTOFF ) 
          {
            overrideAngle = sv.angle;
totalUnboundNucleationCount++;
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
      if (p.useAbsoluteBidirectionalNucBias)
      {
        if (randomGen.randInt(1))
        {
          overrideAngle = p.absoluteNucBiasAngle + randomGen.randNorm(0., p.nucBiasVariance);
        }
        else
        {
          overrideAngle = p.absoluteNucBiasAngle + PI + randomGen.randNorm(0., p.nucBiasVariance);
        }
        while (overrideAngle > 2*PI)
          overrideAngle -= 2*PI;
        while (overrideAngle < 0)
          overrideAngle += 2*PI;
      }
      else
      {
        overrideAngle = handleBiasedNucleationAngle(sv, overrideAngle, preSeeded);
      }
      break;
  }	

  sv.angle = overrideAngle;

  tv = geometry->createTrajectory(sv);
  Microtubule* mt = growing_mts.create(this, tv);
  totalNucleationCount++;
  if (p.useDoubleNucleationSaturation and (not preSeeded))
  {
    occupiedNCs.push_back(systemTime + systemTimeOffset);
  }
  if (p.outputNucPos == "X" or p.outputNucPos == "XY")
  {
    double xPosMeas;
    xPosMeas = sv.x;
    xPosMeas = geometry->xPosGridToDomain(xPosMeas, sv.region->geometryRegionIndex);
    nucleationXpositions.push_back(xPosMeas);
  }
  if (p.outputNucPos == "Y" or p.outputNucPos == "XY")
  {
    double yPosMeas;
    yPosMeas = sv.y;
    yPosMeas = geometry->yPosGridToDomain(yPosMeas, sv.region->geometryRegionIndex);
    nucleationYpositions.push_back(yPosMeas);
  }
  return;
}




void System::determineStochasticEvent()
{

  double nextEventInterval;
  double invTotalRate;
  double typeSelector;
  double logRandom;
  double temp;
  double localCat;
#ifdef BAND_CAT
  localCat = p.catMax;
#else
  localCat = p.kCat;
#endif

  if ((growing_mts.size() == 0) || (p.restrictedPool == 0))
  {
    double constRate;
    double slopeRate;

#ifdef VAR_CAT
    constRate = seedsLeft*p.preSeededRate +p.kNuc*geometry->area + (p.kRes+p.kResExtraMax)*shrinking_mts.size()\
                + p.kSev*totalLength 
                + p.kCat*(growingTipsNormal*(1.+p.kCatAlpha) + p.catastropheMultiplier*growingTipsSpecial)
                + p.kCross*OccupiedIntersectionList.size();
#elif defined (BAND_CAT)
    constRate = seedsLeft*p.preSeededRate +p.kNuc*geometry->area + (p.kRes+p.kResExtraMax)*shrinking_mts.size()\
                + p.kSev*totalLength 
                + p.catMax*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial)
                + p.kCross*OccupiedIntersectionList.size();
#else
    constRate = seedsLeft*p.preSeededRate +p.kNuc*geometry->area + (p.kRes+p.kResExtraMax)*shrinking_mts.size()\
                + p.kSev*totalLength 
                + p.kCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial)
                + p.kCross*OccupiedIntersectionList.size();
#endif
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
  }
  else
  {
    double invLMax = 1.0/(p.poolDensity*geometry->area);
    double alpha = (p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size();
    double beta = p.vPlus*growing_mts.size()*invLMax;
    double L_inf = alpha/beta;
    double baseRate = seedsLeft*p.preSeededRate + p.kNuc*geometry->area + localCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) 
      + p.kRes*shrinking_mts.size() + p.kCross*OccupiedIntersectionList.size();

    double x = p.kSev*(totalLength - L_inf)/(baseRate +  p.kSev*L_inf);
		double y = -beta*log(randomGen.randDblExc())/(baseRate +  p.kSev*L_inf);

		nextEventInterval = (-x+y+ LambertW(x*exp(x-y)))/beta;

		double totalRate = baseRate + p.kSev*(L_inf + (totalLength - L_inf)*exp(-beta*nextEventInterval));
		invTotalRate = 1.0/totalRate;

  }

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
	else if (typeSelector < (p.kNuc*geometry->area + localCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial))*invTotalRate)
		nextStochasticEventType = catastrophe;
	else if (typeSelector < (p.kNuc*geometry->area + localCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + p.kRes*shrinking_mts.size())*invTotalRate)
		nextStochasticEventType = rescue;
	else if (typeSelector < (p.kNuc*geometry->area + localCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + (p.kRes+p.kResExtraMax)*shrinking_mts.size())*invTotalRate)
		nextStochasticEventType = extraRescue;
	else if (typeSelector < (p.kNuc*geometry->area + localCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + (p.kRes+p.kResExtraMax)*shrinking_mts.size() + seedsLeft*p.preSeededRate)*invTotalRate)
		nextStochasticEventType = preSeededNucleation;
	else if (typeSelector < (p.kNuc*geometry->area + localCat*(growingTipsNormal + p.catastropheMultiplier*growingTipsSpecial) + (p.kRes+p.kResExtraMax)*shrinking_mts.size() + seedsLeft*p.preSeededRate + p.kCross*OccupiedIntersectionList.size())*invTotalRate)
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
#ifndef SLEEZY
				exit(-1);
#endif
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

			//flushAndReload(false);		// check out what to do with parameter change events.
			flushAndReload(true);		// check out what to do with parameter change events.

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
        // update double nucleation parameters
        if (p.useDoubleNucleationSaturation)
        {
          initDoubleSatPars();
        }

        // Update velocities and distanceScaleFactors
        Microtubule* mt;
        Segment* seg;
        mt = growing_mts.first();
        while (mt != NULL)
        {
          // always unregister from region (just to be sure)
          mt->plus.trajectory->base.region->unregisterFromRegion(mt->plus.regionTag,t_plus,mt->type);
          mt->minus.trajectory->base.region->unregisterFromRegion(mt->minus.regionTag,t_minus,mt->type);
          mt->plus.velocity = p.vPlus;
          mt->minus.velocity = -p.vTM;
          // reinitialize events and reregister tips on region
          mt->plus.event.reinitialize(mt->plus.event.queue ,mt->plus.velocity);  
          mt->plus.regionTag = mt->plus.trajectory->base.region->registerOnRegion(&(mt->plus),t_plus,mt->type);
          mt->minus.event.reinitialize(mt->minus.event.queue ,mt->minus.velocity); 
          mt->minus.regionTag = mt->minus.trajectory->base.region->registerOnRegion(&(mt->minus),t_minus,mt->type);
          mt->disappearEvent.reinitialize(&(timeQueue), - p.vMin + p.vTM); 
#ifdef DBG_VELOCITY
          //if (randomGen() < 0.1 ) 
          cout << "TEST " << mt->plus.velocity << "\t" << mt->minus.velocity << "\n";
#endif
          mt = mt->next();
        }
        mt = shrinking_mts.first();
        while (mt != NULL)
        {
          // always unregister from region (just to be sure)
          mt->plus.trajectory->base.region->unregisterFromRegion(mt->plus.regionTag,t_plus,mt->type);
          mt->minus.trajectory->base.region->unregisterFromRegion(mt->minus.regionTag,t_minus,mt->type);
          mt->plus.velocity = p.vMin;
          mt->minus.velocity = -p.vTM;
          // reinitialize events and reregister tips on region
          mt->plus.event.reinitialize(mt->plus.event.queue ,mt->plus.velocity);  
          mt->plus.regionTag = mt->plus.trajectory->base.region->registerOnRegion(&(mt->plus),t_plus,mt->type);
          mt->minus.event.reinitialize(mt->minus.event.queue ,mt->minus.velocity);  
          mt->minus.regionTag = mt->minus.trajectory->base.region->registerOnRegion(&(mt->minus),t_minus,mt->type);
          mt->disappearEvent.reinitialize(&(timeQueue), - p.vMin + p.vTM); 
#ifdef DBG_VELOCITY
          //if (randomGen() < 0.1 ) 
          cout << "TEST " << mt->plus.velocity << "\t" << mt->minus.velocity << "\n";
#endif

          mt = mt->next();
        }
        flushAndReload(false);		// RELOAD AGAIN TO FIX QUEUES.




        // Compute pre-seeded nucleation events.
        seedsLeft += static_cast<int> (0.499 + p.preSeededSeedDensity * geometry->area) ;
        // Note that if forbiddenZones==1, a fraction of these events will be 'non-events'


				//// perform further initialization  // Why these 4 lines??
				//mtLengthHistogram.setBins(p.hiresLengthHistogramBins);
				//mtLifetimeHistogram.setBins(p.hiresLifetimeHistogramBins);
				//segAngleLengthHistogram.setBins(p.loresAngleHistogramBins, p.loresLengthHistogramBins);
				//segAngleLifetimeHistogram.setBins(p.loresAngleHistogramBins, p.loresLifetimeHistogramBins);

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
	}    //end of switch statement

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
	is.mirror->ISREF.occupiedListPtr = NULL;
	return;
}

void System::addOccupiedIntersection(IntersectionItr is)
{
	//	cout << "adding OccupiedIntersection\n";
	is->ISREF.occupiedListPtr = OccupiedIntersectionList.create(is);
	is->ISREF.mirror->ISREF.occupiedListPtr = is->ISREF.occupiedListPtr;
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
