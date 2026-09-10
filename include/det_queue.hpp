#ifndef DETERMINISTIC_QUEUE
#define DETERMINISTIC_QUEUE

#include "types.hpp"
#include "parameters.hpp"

class DeterministicQueue
//'Smart' queue object that contains logic to transform to and from the system
// time.
{

  public:

    System* const system; // pointer to containing system
    DeterministicQueue(System*, double (System::*dtFunc)(double), double (System::*tdFunc)(double));

    double currentPos()
    {
        return currentBase; // returns the current 'distance'
    }

    void advanceTime(double); // advances the parameters to the current system time
    void storeTime(int);      // store 'currentBase' in the cache at a given tag position

    double progression(int cachePos) { return currentBase - valueCache[cachePos]; }

    // returns the progress of 'currentBase' relative to a given time tag
    double firstEventTime();  // returns the system time of the first scheduled
                              // event
    DeterministicEvent pop(); // returns and removes the first scheduled event

    bool empty()
    {
        return queue.empty(); // checks whether the queue is empty
    }

    void flush();                             // empties the queue and resets currentBase
    void pushGlobal(double, GlobalEventType); // pushes a global event onto the
                                              // queue at a certain distance
    EventTrackingTag pushDeterministic(double, EventDescriptorIndex);
    // pushes a deterministic (MT) event onto the queue at a certain distance

  private:

    friend class System;
    priority_queue<DeterministicEvent> queue; // the queue itself
    double currentBase;                       // the current 'distance'
    double valueCache[POSITION_CACHE_SIZE];   // cache of previous 'distance' values
    double (System::* const distanceTimeConversionFunction)(double);
    // pointer to a function that converts a distance to a time
    double (System::* const timeDistanceConversionFunction)(double);
    // pointer to a function that converts a time to a distance
};

#endif // DETERMINISTIC_QUEUE
