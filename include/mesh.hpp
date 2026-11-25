#include "types.hpp"
#include "geometry.hpp"
#include "edge.hpp"
#include "surface.hpp"
#include "linalg.hpp"

class Triangle3D
{
  public:

    int element_id;
    int edge_tag;
    int face_tag;
    int poly_intersect_mark;
    int periodic_tag;

    double area;
    double peri_meter;
    double edg_ang_norm;
    double pcat_edge_weight;

    Eigen::Vector3d mid_point;
    Eigen::Vector3d given_normal;

    vector<Edge> side;
    vector<int> vertex_ids;
    vector<int> intersect_edg;
    vector<Eigen::Vector3d> vertex;

    map<int, int> side_map;
    map<int, int> side_rev_map;
    map<int, int> vertex_ids_map;
    map<int, int> side_to_edge;

    Triangle3D()
    {
        element_id = 0;
        edge_tag = 0;
        face_tag = 1;
        poly_intersect_mark = 0;
        periodic_tag = 0;

        area = 0.000001;
        peri_meter = 0.000001;
        edg_ang_norm = 0.000001;
        pcat_edge_weight = 0.000001;
        mid_point << 0.0, 0.0, 0.0;
        given_normal << 0.0, 0.0, 0.0;

        for (int i = 0; i < 3; i++)
        {
            side.push_back(Edge());
            vertex_ids.push_back(0);
            vertex.push_back(Eigen::Vector3d(0, 0, 0));
        }
    }
};

struct EdgeRecord
{
    vector<int> key;
    int cell_id;
    int edge_index;
    EdgeRecord* next;
};

class EdgeOfficer
{
  public:

    int cur;
    double mean_edg_angle;
    vector<int> elements;
    vector<int> edge_sp_store;
    vector<vector<int>> edgeface_store;

    EdgeOfficer(int c)
    {
        cur = c;
        mean_edg_angle = 0.0;
        for (int i = 0; i < 2; i++)
        {
            edgeface_store.push_back(vector<int>());
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
