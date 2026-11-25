#include "types.hpp"
#include "geometry.hpp"
#include "edge.hpp"
#include "surface.hpp"
#include "linalg.hpp"

class Triangle3D
{
  public:

    int elementId;
    int edgeTag;
    int faceTag;
    int polyIntersectMark;
    int periodicTag;

    double area;
    double periMeter;
    double edgAngNorm;
    double pcatEdgeWeight;

    Eigen::Vector3d midPoint;
    Eigen::Vector3d givenNormal;

    vector<Edge> side;
    vector<int> vertexIds;
    vector<int> intersectEdg;
    vector<Eigen::Vector3d> Vertex;

    map<int, int> sideMap;
    map<int, int> sideRevMap;
    map<int, int> vertexIdsMap;
    map<int, int> sideToEdge;

    Triangle3D()
    {
        elementId = 0;
        edgeTag = 0;
        faceTag = 1;
        polyIntersectMark = 0;
        periodicTag = 0;

        area = 0.000001;
        periMeter = 0.000001;
        edgAngNorm = 0.000001;
        pcatEdgeWeight = 0.000001;
        midPoint << 0.0, 0.0, 0.0;
        givenNormal << 0.0, 0.0, 0.0;

        for (int i = 0; i < 3; i++)
        {
            side.push_back(Edge());
            vertexIds.push_back(0);
            Vertex.push_back(Eigen::Vector3d(0, 0, 0));
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
    double meanEdgAngle;
    vector<int> elements;
    vector<int> edgeSpStore;
    vector<vector<int>> edgefaceStore;

    EdgeOfficer(int c)
    {
        cur = c;
        meanEdgAngle = 0.0;
        for (int i = 0; i < 2; i++)
        {
            edgefaceStore.push_back(vector<int>());
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
