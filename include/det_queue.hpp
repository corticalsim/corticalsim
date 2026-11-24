#ifndef DETERMINISTIC_QUEUE
#define DETERMINISTIC_QUEUE

#include "types.hpp"
#include "parameters.hpp"

class DeterministicQueue
//'Smart' queue object that contains logic to transform to and from the system
// time.
{

  public:

    System* const m_system; // pointer to containing system
    DeterministicQueue(System*, double (System::*dt_func)(double), double (System::*td_func)(double));

    double current_pos()
    {
        return m_current_base; // returns the current 'distance'
    }

    void advance_time(double); // advances the parameters to the current system time
    void store_time(int);      // store 'currentBase' in the cache at a given tag position

    double progression(int cache_pos) { return m_current_base - m_value_cache[cache_pos]; }

    // returns the progress of 'currentBase' relative to a given time tag
    double first_event_time();  // returns the system time of the first scheduled
                              // event
    DeterministicEvent pop(); // returns and removes the first scheduled event

    bool empty()
    {
        return m_queue.empty(); // checks whether the queue is empty
    }

    void flush();                             // empties the queue and resets currentBase
    void push_global(double, GlobalEventType); // pushes a global event onto the
                                              // queue at a certain distance
    EventTrackingTag push_deterministic(double, EventDescriptorIndex);
    // pushes a deterministic (MT) event onto the queue at a certain distance

  private:

    friend class System;
    priority_queue<DeterministicEvent> m_queue; // the queue itself
    double m_current_base;                       // the current 'distance'
    double m_value_cache[POSITION_CACHE_SIZE];   // cache of previous 'distance' values
    double (System::* const m_distance_time_conversion_function)(double);
    // pointer to a function that converts a distance to a time
    double (System::* const m_time_distance_conversion_function)(double);
    // pointer to a function that converts a time to a distance
};

#endif // DETERMINISTIC_QUEUE
