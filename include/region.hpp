#ifndef REGION_HPP
#define REGION_HPP

#include "types.hpp"
#include "edge.hpp"
#include "geometry.hpp"
#include "linalg.hpp"

class Region
// virtual base class for a geometry patch with a 2D coordinate system. Within
// a region, region functions are only called by geometry or trajectory member
// functions
{
    friend class Geometry;
    friend class Trajectory;

  public:

    // pointer to geometry
    Geometry* const m_geometry;

    // const RegionType type;
    double m_area;
    double m_total_length;
    // time tag at which the length was last updated
    int m_previous_update_tag;
    int m_region_id;
    double m_z_offset;
    double m_rot_angle;
    double m_peri_meter;
    Eigen::Vector2d m_mid_point;
    vector<Edge> m_side;
    map<int, int> m_side_map;
    map<int, int> m_side_rev_map;
    Eigen::Quaternion<double> m_q;
    vector<Eigen::Vector2d> m_vertices;
    vector<Eigen::Vector3d> m_orientation;
    int m_face_tag;
    int m_poly_intersect_mark;
    vector<int> m_intersect_edg;

    // list of trajectories on the region
    DLList<Trajectory, TRAJECTORY_GRANULARITY> m_trajectories;

    // tip management functions
    list<MTTip*> m_growing_plus_tip_list;
    list<MTTip*> m_shrinking_plus_tip_list;
    list<MTTip*> m_minus_tip_list;
    RegionMTTipTag register_on_region(MTTip*, TipType, MTType);
    void unregister_from_region(RegionMTTipTag, TipType, MTType);

    // updates the length to the current system time tag
    void update_region_length(bool force_update = false);

    Region(Geometry* g, double a):
        m_geometry(g),
        m_area(a),
        m_total_length(0),
        m_previous_update_tag(0),
        m_region_id(0),
        m_z_offset(0.0),
        m_rot_angle(0.0),
        m_peri_meter(0.0)
    {
        m_mid_point << 0.0, 0.0;
        for (int i = 0; i < 3; i++)
        {
            m_vertices.push_back(Eigen::Vector2d(0, 0));
            m_side.push_back(Edge());
        }
    }

    virtual ~Region() { return; };

    double optical_length();

    // vector computation functions (to be implemented in specializations)
    virtual void translate_vector(SurfaceVector&, const double) = 0;
    virtual void get_trajectory_coordinates(SurfaceVector&, double&, vector<PointATedge>&, TrajectoryVector&) = 0;
    virtual double intersection_angle(Trajectory*, Trajectory*) = 0;
    virtual void output_snapshot(ostream&) = 0;
    virtual void output_order_heat_map(ostream&, vector<double>&, vector<Eigen::Vector3d>&) = 0;

  protected:

    // trajectory management functions
    TrajectoryVector insert_trajectory(const SurfaceVector&);
    void remove_trajectory(Trajectory*);
    virtual void make_intersection_list(Trajectory*) = 0;
    virtual SurfaceVector random_surface_vector() = 0;

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

    virtual SurfaceVector random_surface_vector() = 0;
    virtual void get_trajectory_coordinates(SurfaceVector&, double&, vector<PointATedge>&, TrajectoryVector&) = 0;
    void make_intersection_list(Trajectory*);
    void translate_vector(SurfaceVector&, const double);
    double intersection_angle(Trajectory*, Trajectory*);

    void get_order_parameters(OrderParameters&);
    void get_order_parameters_raw_flat(OrderParametersRaw&, vector<Eigen::Vector3d>&, double);
    void get_order_parameters_raw_cylinder(OrderParametersRaw&, double, double, double);

    void output_snapshot(ostream&);
    void output_snapshot_offset(ostream&, double, double);

    void output_order_heat_map(ostream&, vector<double>&, vector<Eigen::Vector3d>&);
};

class Triangle: public Cartesian
{
  public:

    Triangle(double, Geometry*);

    SurfaceVector random_surface_vector();
    void get_trajectory_coordinates(SurfaceVector&, double&, vector<PointATedge>&, TrajectoryVector&);
};

#endif // REGION_HPP
