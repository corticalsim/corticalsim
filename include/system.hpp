#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "randhub.hpp"
#include "types.hpp"
#include "CompactList.hpp"
#include "det_queue.hpp"
#include "event.hpp"
#include "measurement.hpp"
#include "parameters.hpp"
#include "microtubule.hpp"

class System
// Master object defining an interacting mt system
{
  public:

    Geometry* geometry;
    DLList<Microtubule, MICROTUBULE_GRANULARITY> growing_mts;
    DLList<Microtubule, MICROTUBULE_GRANULARITY> shrinking_mts;

#ifdef CROSS_SEV
    CompactList<OccupiedIntersection, OCCUPIED_INTERSECTION_GRANULARITY> OccupiedIntersectionList;
#endif

    vector<int> growingTipsReg;

    Parameters p;

    // random generator
    RandHub<> randomGen;

    // current system time with respect to offset
    double systemTime;

    // current system time offset
    double systemTimeOffset;

    int currentTimeTag;
    double totalLength;
    int memUsage;
    bool stopSignal;

    time_t wallClockStartTime;
    int totalSEventCount;
    int totalValidDEventCount;
    int totalInvalidDEventCount;
    int totalLengthSeveringCount;
    int totalIntersectionSeveringCount;
    int boundaryCrossingCount;
    int totalZipperCount;
    int totalCrossoverCount;
    int totalInducedCatastropheCount;

    int countSegments;
    int countTrajectories;
    int countIntersections;
    int estimateMemoryFootprint();

    double nextStochasticEventTime;
    StochasticEventType nextStochasticEventType;

    // queue object for time-defined events
    DeterministicQueue timeQueue;

    // queue object for vPlus-defined events (regulated by tubulin pool size)
    DeterministicQueue vPlusQueue;

    double identity(double i) { return i; }

    // distance to time conversion for growing tips
    double vPlusToTime(double);

    // time to distance conversion for growing tips
    double timeToVPlus(double);

    double nextStatusEventTime;
    double nextSnapshotEventTime;
    double nextParameterEventTime;
    double nextPPBEventTime;

    // current tag for event descriptor
    EventDescriptorIndex EventDescriptorID;

    // structure that maps EventDescriptorIDs onto EventDescriptor objects
    map<EventDescriptorIndex, EventDescriptor*> EventDescriptorMap;

    // returns a new event descriptor tag
    EventDescriptorIndex registerEventDescriptor(EventDescriptor*);

    // invalidates an event descriptor tag
    void unregisterEventDescriptor(EventDescriptorIndex);

    // returns a pointer to an EventDescriptor through its tag
    EventDescriptor* getEventDescriptor(EventDescriptorIndex);

    // the current event tracking tag
    EventTrackingTag eventID;

    // returns a new event tracking tag
    EventTrackingTag getEventTag() { return eventID++; }

    System(char*);
    ~System();
    bool integrityCheck();
    void flushAndReload(bool = true);
    void advanceTime(double nextTime);

    // updates all positions in the system to the system time
    void updateAll(bool forceUpdate = false);
    void emergencyBreak();

    void run(double, string&, string&);
    void nextEvent(void);
    void handleGlobalEvent(DeterministicEvent&);
    void determineStochasticEvent();
    void handleNucleationEvent();
    void handleSeveringEvent();

#ifdef CROSS_SEV
    void handleSeveringAtCrossEvent();
#endif

    void handleRescueEvent();
    void handleCatastropheEvent();
    void randomPositionOnMicrotubule(double& cutLength, Segment*& cutSeg);

    // returns the result of a collision at a specified angle, and bundle
    // occupancies
    CollisionType collisionResult(double, int, int);

    // returns the probabilities for induced catastrophes and zippering, as a
    // function of the angle
    void collisionProbabilities(double, double&, double&);

    // binomial table. Necessary for bundle collision type: multi-collisions
    double binomialTable[MAXBINOM][MAXBINOM];
    void makeBinomialTable(void);

    double multiPcross(int Npar, int Ncoll, double xSingle, double zSingle);
    double multiPzip(int Npar, int Ncoll, double xSingle, double zSingle);

    ofstream parameterFile;
    ofstream movieFile;
    ofstream heatMapFile;
    ofstream measurementFile;
    ofstream orderDirectorFile;
    deque<Measurement> measurementHistory;
    bool dataSaved;
    void initializeOutput();
    void performMeasurement();
    void writeMeasurementsToFile(int = 0);
    void closeFiles();

#ifdef CROSS_SEV
    // removes pointers from is and its mirror
    void removeOccupiedIntersection(Intersection& is);

    // creates pointers for is and its mirror [occupiedIntersection contains
    // only one pointer: to the one "on top" (to be cut)]
    void addOccupiedIntersection(IntersectionItr is);
#endif

  private:

    // avoid accidental (expensive) copying of Trajectory objects, by declaring
    // private copy constructors without definitions
    System(const System&);
    System& operator=(const System&);
};

#endif // SYSTEM_HPP
