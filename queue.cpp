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

/************************** DeterministicQueue functions ******************/

DeterministicQueue::DeterministicQueue(System* s, double (System::*dtFunc)(double), double (System::*tdFunc)(double)) :
		system(s),
		distanceTimeConversionFunction(dtFunc),
		timeDistanceConversionFunction(tdFunc),
		currentBase(0)
{
// WARNING: Do not use the system pointer in this function, as it has not been initialized
	valueCache[0] = 0;
}


double DeterministicQueue::firstEventTime()
{
	DeterministicEvent event = queue.top();

	#ifdef DBG_ASSERT
	if (event.eventTimeDist - currentBase < -ZERO_CUTOFF)
	{
		cerr << "Queueing error: first event is in the past.\n";
		cerr << "index=" << event.infoIdx << ", global type=" << event.global_type << "\n";
		exit(-1);
	}
	#endif

#ifdef DBG_ACID_TEST
	double diff = (system->*timeDistanceConversionFunction)((system->*distanceTimeConversionFunction)(event.eventTimeDist - currentBase))- (event.eventTimeDist - currentBase);
	if (abs(diff) > 10*ZERO_CUTOFF)
	{
		cerr << "Inversion problem! difference=" << diff << "\n";
	}
#endif

	return system->systemTime + (system->*distanceTimeConversionFunction)(event.eventTimeDist - currentBase);
}

void DeterministicQueue::flush()
{
	while (!queue.empty())
		queue.pop();
	currentBase = 0;	
	valueCache[system->currentTimeTag] = 0;
	return;	
}

DeterministicEvent DeterministicQueue::pop()
{
	DeterministicEvent event = queue.top();
	queue.pop();
	return event;
}

	
void DeterministicQueue::advanceTime(double t)
{
	currentBase += (system->*timeDistanceConversionFunction)(t - system->systemTime);
	return;
}

void DeterministicQueue::storeTime(int tag)
{
	valueCache[tag] = currentBase;
	return;	
}


void DeterministicQueue::pushGlobal(double timedist, GlobalEventType type)
{
	DeterministicEvent event;
	event.infoIdx = -1;
	event.eventTimeDist = timedist;
	event.global_type = type;
	queue.push(event);
	return;	
}

EventTrackingTag DeterministicQueue::pushDeterministic(double timedist, EventDescriptorIndex idx)
{
	DeterministicEvent event;
	event.infoIdx = idx;
	event.eventTimeDist = timedist;
	event.tag = system->getEventTag();
	queue.push(event);
	return event.tag;	
}


/************************ EventDescriptor functions ********************/

EventDescriptor::EventDescriptor(Microtubule* m, DeterministicQueue* q, double velocity) : 
		mt(m), 
		index(0), 
		type(ev_none), 
		tag(-1),
		queue(q)
// on creation, register the object as a valid event descriptor
{
	// WARNING: Do not use the mt pointer in this function, as it has not been initialized
	index = q->system->registerEventDescriptor(this);
	if (abs(velocity) > ZERO_CUTOFF)
		distanceScaleFactor = 1.0/velocity;
	else
	{
		if (velocity <= 0.0)		// inclusive 0: treat v=0 as slow shrinking state
			distanceScaleFactor = -VERY_LARGE;
		else 
			distanceScaleFactor = VERY_LARGE;
	}

	return;
}


EventDescriptor::~EventDescriptor()
// on destruction, remove it again
{
	mt->system->unregisterEventDescriptor(index);
	return;
}

void EventDescriptor::reinitialize(DeterministicQueue* q, double velocity)
{
	queue = q; 
	if (abs(velocity) > ZERO_CUTOFF)
		distanceScaleFactor = 1.0/velocity;
	else
	{
		if (velocity <= 0.0)		// inclusive 0: treat v=0 as slow shrinking state
			distanceScaleFactor = -VERY_LARGE;
		else 
			distanceScaleFactor = VERY_LARGE;
	}
	return;	
}


void EventDescriptor::pushOnQueue(double dist, DeterministicEventType t)
{
	#ifdef DBG_ASSERT
	if (dist*distanceScaleFactor < - 0.0)
	{
		cerr << "DBG/ASSERT: ERROR: creating an event in the past.\n";
		cerr << "system time=" << mt->system->systemTime << ", reference time=" << queue->currentPos() << ", offset=" <<  dist*distanceScaleFactor << ", event type=" << t << ", tag=" << tag << ".\n";
	}
	#endif

	type = t;
	tag = queue->pushDeterministic(queue->currentPos() + dist*distanceScaleFactor, index);
	return;
}
	
void EventDescriptor::clear()
{
	tag = -1;
	type = ev_none;
	return;
}
