#include "types.hpp"
#include "geometry.hpp"
#include "edge.hpp"
#include "surface.hpp"
#include "linalg.hpp"

class Triangle3D
{
  public:

    int m_element_id;
    int m_edge_tag;
    int m_face_tag;
    int m_poly_intersect_mark;
    int m_periodic_tag;

    double m_area;
    double m_peri_meter;
    double m_edg_ang_norm;
    double m_pcat_edge_weight;

    Eigen::Vector3d m_mid_point;
    Eigen::Vector3d m_given_normal;

    vector<Edge> m_side;
    vector<int> m_vertex_ids;
    vector<int> m_intersect_edg;
    vector<Eigen::Vector3d> m_vertex;

    map<int, int> m_side_map;
    map<int, int> m_side_rev_map;
    map<int, int> m_vertex_ids_map;
    map<int, int> m_side_to_edge;

    Triangle3D()
    {
        m_element_id = 0;
        m_edge_tag = 0;
        m_face_tag = 1;
        m_poly_intersect_mark = 0;
        m_periodic_tag = 0;

        m_area = 0.000001;
        m_peri_meter = 0.000001;
        m_edg_ang_norm = 0.000001;
        m_pcat_edge_weight = 0.000001;
        m_mid_point << 0.0, 0.0, 0.0;
        m_given_normal << 0.0, 0.0, 0.0;

        for (int i = 0; i < 3; i++)
        {
            m_side.push_back(Edge());
            m_vertex_ids.push_back(0);
            m_vertex.push_back(Eigen::Vector3d(0, 0, 0));
        }
    }
};

struct EdgeRecord
{
    vector<int> m_key;
    int m_cell_id;
    int m_edge_index;
    EdgeRecord* m_next;
};

class EdgeOfficer
{
  public:

    int m_cur;
    double m_mean_edg_angle;
    vector<int> m_elements;
    vector<int> m_edge_sp_store;
    vector<vector<int>> m_edgeface_store;

    EdgeOfficer(int c)
    {
        m_cur = c;
        m_mean_edg_angle = 0.0;
        for (int i = 0; i < 2; i++)
        {
            m_edgeface_store.push_back(vector<int>());
        }
    }
};

vector<string> split(string, char);

double find_norm(Eigen::VectorXd&);

bool check_eular(int, int, int);

void pickup_shape(Geometry*);

bool check_surface_orientation(vector<Triangle3D*>&);

void image3d_to2_d(vector<Region*>&, vector<Triangle3D*>&, double[]);

void orient_surface(vector<Vertices*>&, vector<Triangle3D*>&, vector<Triangle3D*>&);

void connect_with_globe(vector<Vertices*>&, vector<Triangle3D*>&);

void make_graph(int, vector<ElementList*>&, vector<Triangle3D*>&);

void edge_descriptors(vector<Vertices*>&, vector<ElementList*>&, vector<Triangle3D*>&);

void establish_pbc(vector<Vertices*>&, vector<Triangle3D*>&, vector<Triangle3D*>&, vector<ElementList*>&, double);

void rigid_body_properties(
vector<Vertices*>&, vector<Triangle3D*>&, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&);

void view_graph(vector<Vertices*>&, vector<Triangle3D*>&, Eigen::Vector3d&, string);

bool line_plane_intersect(Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d&);

double area_polygon3_d(vector<Vertices>&, Eigen::Vector3d);

double intersecting_polygon(
vector<Vertices>&, vector<Vertices>&, vector<Region*>&, vector<ElementList>&, Eigen::Vector3d, Eigen::Vector3d, string);
