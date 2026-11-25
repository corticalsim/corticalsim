#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "types.hpp"
#include "surface.hpp"
#include "parameters.hpp"
#include "trajectory.hpp"
#include "linalg.hpp"

class Geometry
// virtual base class for the various types of geometries
{
    friend class System;

  public:

    // const GeometryType type;
    double area;
    double areainMesh;
    System* const system;

    std::vector<ElementList> ppbEdgeList;
    std::vector<Vertices> globalVertex;

    double areaPPB;
    std::vector<Vertices> ppb;

    double patchArea[3];
    std::vector<int> RegionsIndex[3];

    int elementMax;
    Eigen::Vector3d objectCM, objectPA, nucleousPosition;
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

    bool integrity_check();
    double optical_length();
    int trajectory_count();
    void call_translator(SurfaceVector&, int, int);

    virtual void get_order_parameters(OrderParameters&) = 0;
    virtual void output_snapshot(ostream&) = 0;
    virtual void output_order_heat_map(ostream&, vector<double>&, vector<Eigen::Vector3d>&) = 0;

    SurfaceVector random_surface_vector();
    TrajectoryVector create_trajectory(const SurfaceVector&);

    // creates a trajectory at a given surface vector and returns a trajectory
    // vector
    TrajectoryVector create_and_link_trajectory(const SurfaceVector&, Trajectory*, Direction, double, double);

    // only called by Trajectory (friend class)

  private:

    friend TrajectoryVector Trajectory::next_trajectory(Direction);
    virtual TrajectoryVector extend_trajectory(Trajectory*, Direction) = 0;
};

class TriMeshGeometry: public Geometry
{
  public:

    TriMeshGeometry(System*);
    TrajectoryVector extend_trajectory(Trajectory*, Direction);
    void get_order_parameters(OrderParameters&);
    void output_snapshot(ostream&);
    void output_order_heat_map(ostream&, vector<double>&, vector<Eigen::Vector3d>&);
};

#endif // GEOMETRY_HPP
