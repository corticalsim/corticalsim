#ifndef MICROTUBULE_HPP
#define MICROTUBULE_HPP

#include "types.hpp"
#include "event.hpp"
#include "trajectory.hpp"
#include "segment.hpp"
#include "mt_tip.hpp"

class Microtubule: public DLBaseItem<Microtubule>
{
  public:

    System* const m_system;
    MTTip m_plus;
    MTTip m_minus;

    // event descriptor for events that are not directly associated to a single
    // tip (only disappearance)
    EventDescriptor m_disappear_event;
    // type of microtubule (growing, shrinking)
    MTType m_type;
    double m_nucleation_time;

    // time tag at which the positions and lengths were last updated
    int m_previous_update_tag;

    DLList<Segment> m_segments;

    Microtubule(System*, TrajectoryVector, bool = true);
    ~Microtubule();
    bool integrity_check();
    // status functions
    double length();
    void update_length(bool force_update = false);
    // event functions below
    void set_disappear_event();
    void handle_event(const EventDescriptor*);
    void catastrophe();
    void rescue();

    // self-destruct event (plus and minus ends meet)
    void harakiri();
    void wall();
    void collision();
    void zipper(Direction dir);
    void crossover();

    // called when a microtubule retreats across an intersection
    void backtrack(MTTip*);
    void end_of_segment(MTTip*);
    void sever(Segment*, double);
    void sever_at_cross(IntersectionItr is, Segment* cut_seg);
    void split_segment_at_traj_pos(double cut_pos, Segment* cut_seg);

    // translates random position at MT to corresponding position at Segment
    void translate_position_m_t2_segment(double& cut_pos, Segment*& cut_seg);

  private:

    // avoid accidental (disastrous) copying of Trajectory objects, by
    // declaring private copy constructors without definitions
    Microtubule(const Microtubule&);
    Microtubule& operator=(const Microtubule&);
};

#endif // MICROTUBULE_HPP
