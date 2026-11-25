#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include "DLList.hpp"
#include "intersection.hpp"
#include "surface.hpp"

class TrajectoryVector
{
  public:

    double pos;
    Direction dir, Cdir;
    Trajectory* trajectory;

    TrajectoryVector(double p, Direction d, Trajectory* t):
        pos(p),
        dir(d),
        trajectory(t)
    {
    }

    TrajectoryVector() {}

    TrajectoryVector flipped() { return TrajectoryVector(pos, dir == ::FORWARD ? BACKWARD : ::FORWARD, trajectory); }
};

class Trajectory: public DLBaseItem<Trajectory>
// the trajectory is the basic geometrical object. Tips and segments associate
// with a trajectory and move/lie alongside it. Trajectory intersections
// determine the collision points
{
    friend class Region;

  public:

    const SurfaceVector base;
    const double length;

    vector<PointATedge> endPoint;
    TrajectoryVector thisTr;

    // trajectoryVector of the connecting trajectory at the zero end
    TrajectoryVector prevTr;

    // trajectoryVector of the connecting trajectory at the far end
    TrajectoryVector nextTr;

    // cosine of the 3D angle with the previous trajectory (for edge
    // catastrophes)
    double prevTrCosAngle;

    // cosine of the 3D angle with the next trajectory (for edge catastrophes)
    double nextTrCosAngle;

    // regular catastrophe value at the trajectory zero end
    double prevTrpCat;

    // regular catastrophe value at the trajectory far end
    double nextTrpCat;

    multimap<double, Intersection> intersections;

    // sorted list of all intersections
    IntersectionItr wall_end()
    {
        // returns an iterator to the intersection that stands for the far wall
        return intersections.end();
    }

    IntersectionItr wall_begin()
    {
        // returns an iterator to the intersection that stands for the zero
        // wall
        return intersections.begin();
    }

    explicit Trajectory(SurfaceVector, vector<PointATedge>, double);
    ~Trajectory();

    bool integrity_check();

    // returns the connected trajectory in a given direction, and creates it if
    // necessary
    TrajectoryVector next_trajectory(Direction);

    // checks whether the trajectory can safely be removed - and does it
    void conditional_remove();

    // list of pointers to associated segments
    list<Segment*> segments;

    // inserts a segment into the list
    TrjSegmentTag insert_segment(Segment*);

    // removes a segment from the list
    void remove_segment(TrjSegmentTag);

    // list of pointers to associated tips
    list<MTTip*> notificationList;

    // insert a tip for  notifications
    TrjMTTipTag register_for_notifications(MTTip*);

    // removes a tip from the list
    void unregister_for_notifications(TrjMTTipTag);

    // is called whenever an intersection is invalidated by trajectory removal
    void invalidate_intersection(IntersectionItr&);

    // is called when a new intersection is created
    void new_intersection(IntersectionItr&);

    // returns the difference sign of pos1-pos2, or itr1-itr2 if the first,
    // cannot be determined accurately
    int difference_sign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2);

    // total length of trajectory (optical) that is covered by segments
    double covered_length();

    // total length of segments on the trajectory
    double segment_length();

  private:

    // avoid accidental (expensive) copying of Trajectory objects, by declaring
    // private copy constructors without definitions
    Trajectory(const Trajectory&);
    Trajectory& operator=(const Trajectory&);
};

#endif // TRAJECTORY_HPP
