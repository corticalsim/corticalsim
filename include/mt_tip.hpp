#ifndef MT_TIP_HPP
#define MT_TIP_HPP

#include "types.hpp"
#include "event.hpp"

class MTTip
{
  public:

    Microtubule* m_mt;
    Trajectory* m_trajectory;

    Direction m_dir;
    double m_velocity;
    EventDescriptor m_event;
    IntersectionItr m_next_collision;
    TrjMTTipTag m_notification_tag;
    RegionMTTipTag m_region_tag;
    double m_next_event_pos;
    MTTip(Microtubule*, TrajectoryVector, double, DeterministicQueue*, double);
    ~MTTip();
    TipType type();
    double position();
    double other_position();
    Segment& segment();
    void initialize();
    void unlink_from_trajectory();
    void switch_trajectory(Trajectory*, Direction, IntersectionItr, bool = true);
    void locate_intersection();
    void advance_intersection();
    void determine_event();
    void notify_insert(IntersectionItr&);
    void notify_remove(IntersectionItr&);

  private:

    // avoid accidental (expensive) copying of Trajectory objects, by declaring
    // private copy constructors without definitions
    MTTip(const MTTip&);
    MTTip& operator=(const MTTip&);
};

#endif // MT_TIP_HPP
