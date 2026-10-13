#ifndef REGION_HPP
#define REGION_HPP

#include "types.hpp"
#include "edge.hpp"
#include "geometry.hpp"

class Region
// virtual base class for a geometry patch with a 2D coordinate system. Within
// a region, region functions are only called by geometry or trajectory member
// functions
{
    friend class Geometry;
    friend class Trajectory;

  public:

    // pointer to geometry
    Geometry* const geometry;

    // const RegionType type;
    double area;
    double totalLength;
    // time tag at which the length was last updated
    int previousUpdateTag;
    int regionId;
    double zOffset;
    double rotAngle;
    double periMeter;
    Vector2d midPoint;
    vector<Edge> side;
    map<int, int> sideMap;
    map<int, int> sideRevMap;
    Quaternion<double> Q;
    vector<Vector2d> vertices;
    vector<Vector3d> orientation;
    int faceTag;
    int polyIntersectMark;
    vector<int> intersectEdg;

    // list of trajectories on the region
    DLList<Trajectory, TRAJECTORY_GRANULARITY> trajectories;

    // tip management functions
    list<MTTip*> growingPlusTipList;
    list<MTTip*> shrinkingPlusTipList;
    list<MTTip*> minusTipList;
    RegionMTTipTag registerOnRegion(MTTip*, TipType, MTType);
    void unregisterFromRegion(RegionMTTipTag, TipType, MTType);

    // updates the length to the current system time tag
    void updateRegionLength(bool forceUpdate = false);

    Region(Geometry* g, double a):
        geometry(g),
        area(a),
        totalLength(0),
        previousUpdateTag(0),
        regionId(0),
        zOffset(0.0),
        rotAngle(0.0),
        periMeter(0.0)
    {
        midPoint << 0.0, 0.0;
        for (int i = 0; i < 3; i++)
        {
            vertices.push_back(Vector2d(0, 0));
            side.push_back(Edge());
        }
    }

    virtual ~Region() { return; };

    double opticalLength();

    // vector computation functions (to be implemented in specializations)
    virtual void translateVector(SurfaceVector&, const double) = 0;
    virtual void getTrajectoryCoordinates(SurfaceVector&, double&, vector<PointATedge>&, TrajectoryVector&) = 0;
    virtual double intersectionAngle(Trajectory*, Trajectory*) = 0;
    virtual void outputSnapshot(ostream&) = 0;
    virtual void outputOrderHeatMap(ostream&, vector<double>&, vector<Vector3d>&) = 0;

  protected:

    // trajectory management functions
    TrajectoryVector insertTrajectory(const SurfaceVector&);
    void removeTrajectory(Trajectory*);
    virtual void makeIntersectionList(Trajectory*) = 0;
    virtual SurfaceVector randomSurfaceVector() = 0;

  private:

    // avoid accidental (expensive) copying of Region objects, by declaring
    // private copy constructors without definitions
    Region(const Region&);
    Region& operator=(const Region&);
};

class Cartesian: public Region
{
  public:

    Cartesian(Geometry*, double);

    virtual SurfaceVector randomSurfaceVector() = 0;
    virtual void getTrajectoryCoordinates(SurfaceVector&, double&, vector<PointATedge>&, TrajectoryVector&) = 0;
    void makeIntersectionList(Trajectory*);
    void translateVector(SurfaceVector&, const double);
    double intersectionAngle(Trajectory*, Trajectory*);

    void getOrderParameters(OrderParameters&);
    void getOrderParametersRawFlat(OrderParametersRaw&, vector<Vector3d>&, double);
    void getOrderParametersRawCylinder(OrderParametersRaw&, double, double, double);

    void outputSnapshot(ostream&);
    void outputSnapshotOffset(ostream&, double, double);

    void outputOrderHeatMap(ostream&, vector<double>&, vector<Vector3d>&);
};

class Triangle: public Cartesian
{
  public:

    Triangle(double, Geometry*);

    SurfaceVector randomSurfaceVector();
    void getTrajectoryCoordinates(SurfaceVector&, double&, vector<PointATedge>&, TrajectoryVector&);
};

#endif // REGION_HPP
