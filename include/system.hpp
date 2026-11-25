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
#include "mt_tip.hpp"

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
    int estimate_memory_footprint();

    double nextStochasticEventTime;
    StochasticEventType nextStochasticEventType;

    // queue object for time-defined events
    DeterministicQueue timeQueue;

    // queue object for vPlus-defined events (regulated by tubulin pool size)
    DeterministicQueue vPlusQueue;

    double identity(double i) { return i; }

    // distance to time conversion for growing tips
    double v_plus_to_time(double);

    // time to distance conversion for growing tips
    double time_to_v_plus(double);

    double nextStatusEventTime;
    double nextSnapshotEventTime;
    double nextParameterEventTime;
    double nextPPBEventTime;

    // current tag for event descriptor
    EventDescriptorIndex EventDescriptorID;

    // structure that maps EventDescriptorIDs onto EventDescriptor objects
    map<EventDescriptorIndex, EventDescriptor*> EventDescriptorMap;

    // returns a new event descriptor tag
    EventDescriptorIndex register_event_descriptor(EventDescriptor*);

    // invalidates an event descriptor tag
    void unregister_event_descriptor(EventDescriptorIndex);

    // returns a pointer to an EventDescriptor through its tag
    EventDescriptor* get_event_descriptor(EventDescriptorIndex);

    // the current event tracking tag
    EventTrackingTag eventID;

    // returns a new event tracking tag
    EventTrackingTag get_event_tag() { return eventID++; }

    System(char*);
    ~System();
    bool integrity_check();
    void flush_and_reload(bool = true);
    void advance_time(double next_time);

    // updates all positions in the system to the system time
    void update_all(bool force_update = false);
    void emergency_break();

    void run(double, string&, string&);
    void next_event(void);
    void handle_global_event(DeterministicEvent&);
    void determine_stochastic_event();
    void handle_nucleation_event();
    void handle_severing_event();

#ifdef CROSS_SEV
    void handle_severing_at_cross_event();
#endif

    void handle_rescue_event();
    void handle_catastrophe_event();
    void random_position_on_microtubule(double& cut_length, Segment*& cut_seg);

    // returns the result of a collision at a specified angle, and bundle
    // occupancies
    CollisionType collision_result(double, int, int);

    // returns the probabilities for induced catastrophes and zippering, as a
    // function of the angle
    void collision_probabilities(double, double&, double&);

    // binomial table. Necessary for bundle collision type: multi-collisions
    double binomialTable[MAXBINOM][MAXBINOM];
    void make_binomial_table(void);

    double multi_pcross(int npar, int ncoll, double x_single, double z_single);
    double multi_pzip(int npar, int ncoll, double x_single, double z_single);

    ofstream parameterFile;
    ofstream movieFile;
    ofstream heatMapFile;
    ofstream measurementFile;
    ofstream orderDirectorFile;
    deque<Measurement> measurementHistory;
    bool dataSaved;
    void initialize_output();
    void perform_measurement();
    void write_measurements_to_file(int = 0);
    void close_files();

#ifdef CROSS_SEV
    // removes pointers from is and its mirror
    void remove_occupied_intersection(Intersection& is);

    // creates pointers for is and its mirror [occupiedIntersection contains
    // only one pointer: to the one "on top" (to be cut)]
    void add_occupied_intersection(IntersectionItr is);
#endif

  private:

    // avoid accidental (expensive) copying of Trajectory objects, by declaring
    // private copy constructors without definitions
    System(const System&);
    System& operator=(const System&);
};

#endif // SYSTEM_HPP
