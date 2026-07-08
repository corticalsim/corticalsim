#ifndef SEGMENT_HPP
#define SEGMENT_HPP

class Segment: public DLBaseItem<Segment>
{
  public:

    Microtubule* mt;
    Trajectory* const trajectory;
    TrjSegmentTag trajectoryTag;
    double nucleationTime;

    // -+ direction with respect to trajectory orientation
    Direction dir;
    double start;
    double end;
    IntersectionItr startItr;
    IntersectionItr endItr;

    // constructs a segment as part of a microtubule at a specified vector
    // location
    Segment(Microtubule*, TrajectoryVector&);

    ~Segment();

    double length() { return abs(end - start); }

    bool isLastInMT();
    bool isFirstInMT();
    bool crossesIntersection(IntersectionItr& is);
};

#endif // SEGMENT_HPP
