#include "det_queue.hpp"
#include "system.hpp"

DeterministicQueue::DeterministicQueue(System* s, double (System::*dt_func)(double), double (System::*td_func)(double)):
    m_system(s),
    m_current_base(0),
    m_distance_time_conversion_function(dt_func),
    m_time_distance_conversion_function(td_func)
{
    // do not use the system pointer here, as it has not been initialized
    m_value_cache[0] = 0;
}

double DeterministicQueue::first_event_time()
{
    DeterministicEvent event = m_queue.top();

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

    return m_system->m_system_time + (m_system->*m_distance_time_conversion_function)(event.m_event_time_dist - m_current_base);
}

void DeterministicQueue::flush()
{
    // flush queue
    while (!m_queue.empty())
    {
        m_queue.pop();
    }

    m_current_base = 0;
    m_value_cache[m_system->m_current_time_tag] = 0;

    return;
}

DeterministicEvent DeterministicQueue::pop()
{
    // get a deterministic event from top of the queue
    DeterministicEvent event = m_queue.top();
    m_queue.pop();

    return event;
}

void DeterministicQueue::advance_time(double t)
{
    // get advanced in time
    m_current_base += (m_system->*m_time_distance_conversion_function)(t - m_system->m_system_time);
    return;
}

void DeterministicQueue::store_time(int tag)
{
    // store the current base time
    m_value_cache[tag] = m_current_base;
    return;
}

void DeterministicQueue::push_global(double timedist, GlobalEventType type)
{
    DeterministicEvent event;
    event.m_info_idx = -1;
    event.m_event_time_dist = timedist;
    event.m_global_type = type;

    // store a new global event into the queue
    m_queue.push(event);
    return;
}

EventTrackingTag DeterministicQueue::push_deterministic(double timedist, EventDescriptorIndex idx)
{
    DeterministicEvent event;
    event.m_info_idx = idx;
    event.m_event_time_dist = timedist;
    event.m_tag = m_system->get_event_tag();

    // store a new deterministic event into the queue
    m_queue.push(event);
    return event.m_tag;
}

EventDescriptor::EventDescriptor(Microtubule* m, DeterministicQueue* q, double velocity):
    m_index(0),
    m_mt(m),
    m_type(EV_NONE),
    m_tag(-1),
    m_queue(q)
{
    // after creation, register the object as a valid event descriptor
    m_index = q->m_system->register_event_descriptor(this);

    if (abs(velocity) > ZERO_CUTOFF)
    {
        m_distance_scale_factor = 1.0 / velocity;
    }
    else
    {
        if (velocity == 0.0)
        {
            m_distance_scale_factor = -VERY_LARGE;
        }
        else if (velocity >= 0.0)
        {
            m_distance_scale_factor = VERY_LARGE;
        }
        else
        {
            m_distance_scale_factor = -VERY_LARGE;
        }
    }

    return;
}

EventDescriptor::~EventDescriptor()
{
    // unregister the event descriptor
    m_mt->m_system->unregister_event_descriptor(m_index);
    return;
}

void EventDescriptor::reinitialize(DeterministicQueue* q, double velocity)
{
    // re-initialize the event descriptor
    m_queue = q;
    if (abs(velocity) < 1. / ZERO_CUTOFF)
    {
        m_distance_scale_factor = 1.0 / velocity;
    }
    else
    {
        if (velocity >= 0.0)
        {
            m_distance_scale_factor = VERY_LARGE;
        }
        else
        {
            m_distance_scale_factor = -VERY_LARGE;
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

    m_type = t;

    // store deterministic event into the queue
    m_tag = m_queue->push_deterministic(m_queue->current_pos() + dist * m_distance_scale_factor, m_index);
    return;
}

void EventDescriptor::clear()
{
    // clear event descriptor
    m_tag = -1;
    m_type = EV_NONE;
    return;
}
