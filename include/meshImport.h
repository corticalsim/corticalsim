#include "corticalSimReal.h"

using namespace Eigen;
using namespace std;

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
	double  pcatEdgeWeight;

        Vector3d midPoint;
        Vector3d givenNormal;

        vector<Edge> side;
        vector<int> vertexIds;
        vector<int> intersectEdg;
        vector<Vector3d> Vertex;

        map<int,int> sideMap;
        map<int,int> sideRevMap;
        map<int,int> vertexIdsMap;
        map<int,int> sideToEdge;

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
            midPoint << 0.0,0.0,0.0;
            givenNormal << 0.0,0.0,0.0;

            for(int i=0;i<3;i++)
            {
                side.push_back(Edge());
                vertexIds.push_back(0);
                Vertex.push_back(Vector3d(0,0,0));
            }
        }
};

struct edgeRecord {
	vector<int> key;
	int cell_id;
	int edge_index;
	edgeRecord* next;
};

class edgeOfficer
{
	public:
	int cur;
        double meanEdgAngle;
	vector<int> elements;
	vector<int> edgeSpStore;
	vector<vector<int> > edgefaceStore;

	edgeOfficer(int c)
	{
              cur = c;
	      meanEdgAngle = 0.0;
	      for(int i=0;i<2;i++)
	      edgefaceStore.push_back(vector<int>());
	}
};

vector<string> split(string, char);

double findNorm(VectorXd &);

bool checkEular(int,int,int);

void pickup_shape(Geometry*);

bool checkSurfaceOrientation(vector<Triangle3D*>&);

void Image3dTo2D(vector<Region*> &,vector<Triangle3D*> &, double []);

void orientSurface(vector<Vertics*>&,vector<Triangle3D*>&,vector<Triangle3D*>&);

void connectWithGlobe(vector<Vertics*>&,vector<Triangle3D*>&);

void makeGraph(int ,vector<elementList*>&,vector<Triangle3D*>&);

void edgeDescriptors(vector<Vertics*>&,vector<elementList*> &,vector<Triangle3D*>&);

void establishPBC(vector<Vertics*>&,vector<Triangle3D*>&,vector<Triangle3D*>&,vector<elementList*> &,double);

void rigidBodyProperties(vector<Vertics*>&,vector<Triangle3D*>&,Vector3d &,Vector3d &,Vector3d &,Vector3d &);

void viewGraph(vector<Vertics*>&,vector<Triangle3D*>&,Vector3d &,string);

bool linePlaneIntersect(Vector3d, Vector3d, Vector3d, Vector3d, Vector3d&);

double areaPolygon3D(vector<Vertics>&,Vector3d);

double intersectingPolygon(vector<Vertics>&,vector<Vertics>&,vector<Region*>&,vector<elementList> &,Vector3d,Vector3d,string);



