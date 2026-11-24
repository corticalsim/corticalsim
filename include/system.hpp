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

    Geometry* m_geometry;
    DLList<Microtubule, MICROTUBULE_GRANULARITY> m_growing_mts;
    DLList<Microtubule, MICROTUBULE_GRANULARITY> m_shrinking_mts;

#ifdef CROSS_SEV
    CompactList<OccupiedIntersection, OCCUPIED_INTERSECTION_GRANULARITY> m_occupied_intersection_list;
#endif

    vector<int> m_growing_tips_reg;

    Parameters m_p;

    // random generator
    RandHub<> m_random_gen;

    // current system time with respect to offset
    double m_system_time;

    // current system time offset
    double m_system_time_offset;

    int m_current_time_tag;
    double m_total_length;
    int m_mem_usage;
    bool m_stop_signal;

    time_t m_wall_clock_start_time;
    int m_total_s_event_count;
    int m_total_valid_d_event_count;
    int m_total_invalid_d_event_count;
    int m_total_length_severing_count;
    int m_total_intersection_severing_count;
    int m_boundary_crossing_count;
    int m_total_zipper_count;
    int m_total_crossover_count;
    int m_total_induced_catastrophe_count;

    int m_count_segments;
    int m_count_trajectories;
    int m_count_intersections;
    int estimate_memory_footprint();

    double m_next_stochastic_event_time;
    StochasticEventType m_next_stochastic_event_type;

    // queue object for time-defined events
    DeterministicQueue m_time_queue;

    // queue object for vPlus-defined events (regulated by tubulin pool size)
    DeterministicQueue m_v_plus_queue;

    double identity(double i) { return i; }

    // distance to time conversion for growing tips
    double v_plus_to_time(double);

    // time to distance conversion for growing tips
    double time_to_v_plus(double);

    double m_next_status_event_time;
    double m_next_snapshot_event_time;
    double m_next_parameter_event_time;
    double m_next_ppb_event_time;

    // current tag for event descriptor
    EventDescriptorIndex m_event_descriptor_id;

    // structure that maps EventDescriptorIDs onto EventDescriptor objects
    map<EventDescriptorIndex, EventDescriptor*> m_event_descriptor_map;

    // returns a new event descriptor tag
    EventDescriptorIndex register_event_descriptor(EventDescriptor*);

    // invalidates an event descriptor tag
    void unregister_event_descriptor(EventDescriptorIndex);

    // returns a pointer to an EventDescriptor through its tag
    EventDescriptor* get_event_descriptor(EventDescriptorIndex);

    // the current event tracking tag
    EventTrackingTag m_event_id;

    // returns a new event tracking tag
    EventTrackingTag get_event_tag() { return m_event_id++; }

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
    double m_binomial_table[MAXBINOM][MAXBINOM];
    void make_binomial_table(void);

    double multi_pcross(int npar, int ncoll, double x_single, double z_single);
    double multi_pzip(int npar, int ncoll, double x_single, double z_single);

    ofstream m_parameter_file;
    ofstream m_movie_file;
    ofstream m_heat_map_file;
    ofstream m_measurement_file;
    ofstream m_order_director_file;
    deque<Measurement> m_measurement_history;
    bool m_data_saved;
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
