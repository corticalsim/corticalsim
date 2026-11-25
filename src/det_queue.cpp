#include "det_queue.hpp"
#include "system.hpp"

DeterministicQueue::DeterministicQueue(System* s, double (System::*dt_func)(double), double (System::*td_func)(double)):
    system(s),
    currentBase(0),
    distanceTimeConversionFunction(dt_func),
    timeDistanceConversionFunction(td_func)
{
    // do not use the system pointer here, as it has not been initialized
    valueCache[0] = 0;
}

double DeterministicQueue::first_event_time()
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
    double diff = (system->*timeDistanceConversionFunction)(
                  (system->*distanceTimeConversionFunction)(event.eventTimeDist - currentBase))
                  - (event.eventTimeDist - currentBase);
    if (abs(diff) > 10 * ZERO_CUTOFF)
    {
        cerr << "Inversion problem! difference=" << diff << "\n";
    }
#endif

    return system->systemTime + (system->*distanceTimeConversionFunction)(event.eventTimeDist - currentBase);
}

void DeterministicQueue::flush()
{
    // flush queue
    while (!queue.empty())
    {
        queue.pop();
    }

    currentBase = 0;
    valueCache[system->currentTimeTag] = 0;

    return;
}

DeterministicEvent DeterministicQueue::pop()
{
    // get a deterministic event from top of the queue
    DeterministicEvent event = queue.top();
    queue.pop();

    return event;
}

void DeterministicQueue::advance_time(double t)
{
    // get advanced in time
    currentBase += (system->*timeDistanceConversionFunction)(t - system->systemTime);
    return;
}

void DeterministicQueue::store_time(int tag)
{
    // store the current base time
    valueCache[tag] = currentBase;
    return;
}

void DeterministicQueue::push_global(double timedist, GlobalEventType type)
{
    DeterministicEvent event;
    event.infoIdx = -1;
    event.eventTimeDist = timedist;
    event.global_type = type;

    // store a new global event into the queue
    queue.push(event);
    return;
}

EventTrackingTag DeterministicQueue::push_deterministic(double timedist, EventDescriptorIndex idx)
{
    DeterministicEvent event;
    event.infoIdx = idx;
    event.eventTimeDist = timedist;
    event.tag = system->get_event_tag();

    // store a new deterministic event into the queue
    queue.push(event);
    return event.tag;
}

EventDescriptor::EventDescriptor(Microtubule* m, DeterministicQueue* q, double velocity):
    index(0),
    mt(m),
    type(EV_NONE),
    tag(-1),
    queue(q)
{
    // after creation, register the object as a valid event descriptor
    index = q->system->register_event_descriptor(this);

    if (abs(velocity) > ZERO_CUTOFF)
    {
        distanceScaleFactor = 1.0 / velocity;
    }
    else
    {
        if (velocity == 0.0)
        {
            distanceScaleFactor = -VERY_LARGE;
        }
        else if (velocity >= 0.0)
        {
            distanceScaleFactor = VERY_LARGE;
        }
        else
        {
            distanceScaleFactor = -VERY_LARGE;
        }
    }

    return;
}

EventDescriptor::~EventDescriptor()
{
    // unregister the event descriptor
    mt->system->unregister_event_descriptor(index);
    return;
}

void EventDescriptor::reinitialize(DeterministicQueue* q, double velocity)
{
    // re-initialize the event descriptor
    queue = q;
    if (abs(velocity) < 1. / ZERO_CUTOFF)
    {
        distanceScaleFactor = 1.0 / velocity;
    }
    else
    {
        if (velocity >= 0.0)
        {
            distanceScaleFactor = VERY_LARGE;
        }
        else
        {
            distanceScaleFactor = -VERY_LARGE;
        }
    }

    return;
}

void EventDescriptor::push_on_queue(double dist, DeterministicEventType t)
{
#ifdef DBG_ASSERT
    if (dist * distanceScaleFactor < -0.0)
    {
        cerr << "DBG/ASSERT: ERROR: creating an event in the past.\n";
        cerr << "system time=" << mt->system->systemTime << ", reference time=" << queue->currentPos()
             << ", offset=" << dist * distanceScaleFactor << ", event type=" << t << ", tag=" << tag << ".\n";
    }
#endif

    type = t;

    // store deterministic event into the queue
    tag = queue->push_deterministic(queue->current_pos() + dist * distanceScaleFactor, index);
    return;
}

void EventDescriptor::clear()
{
    // clear event descriptor
    tag = -1;
    type = EV_NONE;
    return;
}
