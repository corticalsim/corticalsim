#ifndef MICROTUBULE_HPP
#define MICROTUBULE_HPP

#include "types.hpp"
#include "event.hpp"
#include "det_queue.hpp"
#include "trajectory.hpp"
#include "segment.hpp"

class MTTip
{
  public:

    Microtubule* mt;
    Trajectory* trajectory;

    Direction dir;
    double velocity;
    EventDescriptor event;
    IntersectionItr nextCollision;
    TrjMTTipTag notificationTag;
    RegionMTTipTag regionTag;
    double nextEventPos;
    MTTip(Microtubule*, TrajectoryVector, double, DeterministicQueue*, double);
    ~MTTip();
    TipType type();
    double position();
    double otherPosition();
    Segment& segment();
    void initialize();
    void unlinkFromTrajectory();
    void switchTrajectory(Trajectory*, Direction, IntersectionItr, bool = true);
    void locateIntersection();
    void advanceIntersection();
    void determineEvent();
    void notifyInsert(IntersectionItr&);
    void notifyRemove(IntersectionItr&);

  private:

    // avoid accidental (expensive) copying of Trajectory objects, by declaring
    // private copy constructors without definitions
    MTTip(const MTTip&);
    MTTip& operator=(const MTTip&);
};

class Microtubule: public DLBaseItem<Microtubule>
{
  public:

    System* const system;
    MTTip plus;
    MTTip minus;

    // event descriptor for events that are not directly associated to a single
    // tip (only disappearance)
    EventDescriptor disappearEvent;
    // type of microtubule (growing, shrinking)
    MTType type;
    double nucleationTime;

    // time tag at which the positions and lengths were last updated
    int previousUpdateTag;

    DLList<Segment> segments;

    Microtubule(System*, TrajectoryVector, bool = true);
    ~Microtubule();
    bool integrityCheck();
    // status functions
    double length();
    void updateLength(bool forceUpdate = false);
    // event functions below
    void setDisappearEvent();
    void handleEvent(const EventDescriptor*);
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
    void endOfSegment(MTTip*);
    void sever(Segment*, double);
    void severAtCross(IntersectionItr is, Segment* cutSeg);
    void splitSegmentAtTrajPos(double cutPos, Segment* cutSeg);

    // translates random position at MT to corresponding position at Segment
    void translatePositionMT2Segment(double& cutPos, Segment*& cutSeg);

  private:

    // avoid accidental (disastrous) copying of Trajectory objects, by
    // declaring private copy constructors without definitions
    Microtubule(const Microtubule&);
    Microtubule& operator=(const Microtubule&);
};

#endif // MICROTUBULE_HPP
