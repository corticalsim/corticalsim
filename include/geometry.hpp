#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "types.hpp"
#include "surface.hpp"
#include "parameters.hpp"
#include "trajectory.hpp"

class Geometry
// virtual base class for the various types of geometries
{
    friend class System;

  public:

    // const GeometryType type;
    double area;
    double areainMesh;
    System* const system;

    std::vector<elementList> ppbEdgeList;
    std::vector<Vertices> globalVertex;

    double areaPPB;
    std::vector<Vertices> ppb;

    double patchArea[3];
    std::vector<int> RegionsIndex[3];

    int elementMax;
    Vector3d objectCM, objectPA, nucleousPosition;
    std::vector<Region*> regions;

    Geometry(System* s, double a):
        area(a),
        system(s)
    {
        for (int i = 0; i < 3; i++)
        {
            patchArea[i] = 0.0;
        }
    };

    // virtual destructor implies virtual class
    virtual ~Geometry() { return; };

    bool integrityCheck();
    double opticalLength();
    int trajectoryCount();
    void callTranslator(SurfaceVector&, int, int);

    virtual void getOrderParameters(OrderParameters&) = 0;
    virtual void outputSnapshot(ostream&) = 0;
    virtual void outputOrderHeatMap(ostream&, vector<double>&, vector<Vector3d>&) = 0;

    SurfaceVector randomSurfaceVector();
    TrajectoryVector createTrajectory(const SurfaceVector&);

    // creates a trajectory at a given surface vector and returns a trajectory
    // vector
    TrajectoryVector createAndLinkTrajectory(const SurfaceVector&, Trajectory*, Direction, double, double);

    // only called by Trajectory (friend class)

  private:

    friend TrajectoryVector Trajectory::nextTrajectory(Direction);
    virtual TrajectoryVector extendTrajectory(Trajectory*, Direction) = 0;
};

class TriMeshGeometry: public Geometry
{
  public:

    TriMeshGeometry(System*);
    TrajectoryVector extendTrajectory(Trajectory*, Direction);
    void getOrderParameters(OrderParameters&);
    void outputSnapshot(ostream&);
    void outputOrderHeatMap(ostream&, vector<double>&, vector<Vector3d>&);
};

#endif // GEOMETRY_HPP
