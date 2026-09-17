#ifndef EVENT_HPP
#define EVENT_HPP

#include "types.hpp"
#include "det_queue.hpp"

class EventDescriptor
{
    // The EventDescriptor object is the microtubule-side interface to the
    // event queues. It stores information about the upcoming event, and is
    // used as an entry-point for the event handler when this event is
    // executed. Every microtubule has three of these objects: one on every tip
    // and one to keep track of 'ev_disappear' events, when the two tips
    // annihilate each other.

  public:

    // index of this object in the EventDescriptorMap (necessary for removal)
    EventDescriptorIndex index;
    Microtubule* const mt;

    // type of queued event
    DeterministicEventType type;

    // associated tracking tag (for validity testing)
    EventTrackingTag tag;

    // inverse velocity (with respect to the event clock defined in the queue)
    double distanceScaleFactor;
    DeterministicQueue* queue;
    EventDescriptor(Microtubule*, DeterministicQueue*, double);
    ~EventDescriptor();

    // resets the queue and velocity of the event timer
    void reinitialize(DeterministicQueue*, double);

    // pushes a deterministic event at a certain distance onto the queue
    void pushOnQueue(double, DeterministicEventType);

    // invalidates the queued event
    void clear();
};

class DeterministicEvent
{
  public:

    double eventTimeDist;

    // -1 indicates a global event
    EventDescriptorIndex infoIdx;

    union
    {
        EventTrackingTag tag;
        GlobalEventType global_type;
    };

    // define the < operator for automatic sorting of events in the event queue
    // (nearest first)
    bool operator<(const DeterministicEvent& ev2) const { return eventTimeDist > ev2.eventTimeDist; }
};

#endif // EVENT_HPP
