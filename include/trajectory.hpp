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

    TrajectoryVector flipped() { return TrajectoryVector(pos, dir == ::forward ? backward : ::forward, trajectory); }
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
    IntersectionItr wallEnd()
    {
        // returns an iterator to the intersection that stands for the far wall
        return intersections.end();
    }

    IntersectionItr wallBegin()
    {
        // returns an iterator to the intersection that stands for the zero
        // wall
        return intersections.begin();
    }

    explicit Trajectory(SurfaceVector, vector<PointATedge>, double);
    ~Trajectory();

    bool integrityCheck();

    // returns the connected trajectory in a given direction, and creates it if
    // necessary
    TrajectoryVector nextTrajectory(Direction);

    // checks whether the trajectory can safely be removed - and does it
    void conditionalRemove();

    // list of pointers to associated segments
    list<Segment*> segments;

    // inserts a segment into the list
    TrjSegmentTag insertSegment(Segment*);

    // removes a segment from the list
    void removeSegment(TrjSegmentTag);

    // list of pointers to associated tips
    list<MTTip*> notificationList;

    // insert a tip for  notifications
    TrjMTTipTag registerForNotifications(MTTip*);

    // removes a tip from the list
    void unregisterForNotifications(TrjMTTipTag);

    // is called whenever an intersection is invalidated by trajectory removal
    void invalidateIntersection(IntersectionItr&);

    // is called when a new intersection is created
    void newIntersection(IntersectionItr&);

    // returns the difference sign of pos1-pos2, or itr1-itr2 if the first,
    // cannot be determined accurately
    int differenceSign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2);

    // total length of trajectory (optical) that is covered by segments
    double coveredLength();

    // total length of segments on the trajectory
    double segmentLength();

  private:

    // avoid accidental (expensive) copying of Trajectory objects, by declaring
    // private copy constructors without definitions
    Trajectory(const Trajectory&);
    Trajectory& operator=(const Trajectory&);
};

#endif // TRAJECTORY_HPP
