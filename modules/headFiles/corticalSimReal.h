#ifndef CORTICALSIM_H_
#define CORTICALSIM_H_

/*******************************************************************************************************************************************************/
//For details of simulation implementation:
//Chakrabortty B, Blilou I, Scheres B, Mulder BM (2018) A computational framework for cortical microtubule dynamics in realistically shaped plant cells.
//PLoS Comput Biol 14(2): e1005959. https://doi.org/10.1371/journal.pcbi.1005959

/***********************************************************************************************/
// This is a framework to simulate cortical microtubule dynamics in arbitrary Shaped plant cells.
// The program is developed by Bandan Chakrabortty (07/03/2013) and some basic implementations
// are used from (cortSim:1.20) developed by Simon Tindemans (01/04/2008)

/*********************************************************************************************/
/*
 * IMPORTANT: pitfalls to watch out for when extending the code
 * * The microtubule is only assumed to grow in length at the plus end. This is not an MTTip property, but a
 *   property of the Microtubule functions. Two-sided growth can be enabled, but needs plus/minus distinctions
 *   for the events. Also, tips with velocity zero (usually non-treadmilling minus ends) are
 *   assumed to be in a 'shrinking' state.
 * * The program assumes a stable sort order for the intersection lists. This is true in practice (in all common
 *   implementations), but not mandated by the standard.
 * * If the tubulin pool size is decreased on the fly, care should be taken that the actual density doesn't exceed the
 *   pool size. In this situation, the behaviour is undefined.
 * * In microtubule.cpp and system.cpp, compiler warnings are (should be) issued for using the 'this' pointer within
 *   the initializer list, because the 'this' pointer cannot be used until the initialization has completed. The
 *   current use is ok, because the pointer is only used to store the address of the parent object. NOTE: when making
 *   changes to the code, take care NOT to use this pointer within the constructors of other objects!
 */

#define PROGRAM_VERSION "BC.2018"
#define CROSS_SEV
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <boost/pool/pool_alloc.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "Eigen/Dense"
#include "MersenneTwister.h"
#include "eig3.h"
#include "DLList.h"
#include <float.h>
#include "CompactList.h"
#include <sys/stat.h>

using namespace Eigen;
using namespace std;

// numeric parameters
const double PI = 3.141592653589793;
const double ZERO_CUTOFF = 1000000*numeric_limits<double>::epsilon();   // approx 10E-10;//
const double VERY_LARGE = 10E100;
const int MAXBINOM = 100;

// memory management
const int MICROTUBULE_GRANULARITY = 256;
const int SEGMENT_GRANULARITY = 1024;
const int TRAJECTORY_GRANULARITY = 512;
const int INTERSECTION_GRANULARITY = 4096;
const int EVENT_GRANULARITY = 1024;
const int OCCUPIED_INTERSECTION_GRANULARITY = 256;

// maximum size of measurement cache (increase to minimize disk access)
const int MAX_HISTORY_SIZE = 256;

// minimum size of measurement cache (minimum may be needed for trending)
const int MIN_HISTORY_SIZE = 64;

const int POSITION_CACHE_SIZE = 32768;
const int CLOCK_POLLING_INTERVAL = 10000;
const int MEMORY_POLLING_INTERVAL = 200;

// this is related to the total number of events
const int QUEUE_FLUSH_INTERVAL = 1000000;

// used to create bins for biased nucleation
const int NUCLEATION_DISCRETIZATION_STEPS = 512;

// new data type definition
typedef enum {forward=1, backward=-1} Direction;
extern string DirectionTypeText[];
typedef enum {ev_none, ev_wall=1, ev_collision, ev_backtrack, ev_end_of_segment, ev_disappear} DeterministicEventType;
typedef enum {measure=1, snapshot, stop, status,parameter_change, signalPPB} GlobalEventType;
typedef enum {catastrophe=1, rescue, katanin, severingAtCross, nucleation} StochasticEventType;
typedef enum {ct_zipper=1, ct_crossover, ct_inducedCatastrophe} CollisionType;
typedef enum {nuc_isotropic, nuc_COUNT_LAST} NucleationType;
extern string NucleationTypeText[];
typedef enum {int_zipFirst,int_catFirst,int_COUNT_LAST} InteractionType;
extern string InteractionTypeText[];
typedef enum {bdl_simple,  bdl_sticky, bdl_noZip, bdl_multiCollision, bdl_Ncollision, bdl_COUNT_LAST} BundleType;
extern string BundleTypeText[];
typedef enum {t_minus, t_plus} TipType;
typedef enum {mt_growing, mt_shrinking} MTType;

class System;
class Microtubule;
class Segment;
class MTTip;
class DeterministicEvent;
class Parameters;
class Geometry;
class Region;
class Trajectory;
class Intersection;

#ifdef CROSS_SEV
class OccupiedIntersection;
#endif

typedef multimap<double, Intersection, std::less<double>, boost::fast_pool_allocator<std::pair<double, Intersection> > >::iterator IntersectionItr;
typedef int EventDescriptorIndex;
typedef int EventTrackingTag;
typedef list<Segment*>::iterator TrjSegmentTag;
typedef list<MTTip*>::iterator TrjMTTipTag;
typedef list<MTTip*>::iterator RegionMTTipTag;

class Vertics
{
public:
    double x;
    double y;
    double z;
    Vertics(double x1,double y1,double z1)
    {
        x = x1;
        y = y1;
        z = z1;
    }
};

struct SurfaceVector
{
    double x;
    double y;
    double z;
    double angle;
    double tvPos;

    Region* region;

    SurfaceVector()
    {
        x = 0.0;
        y = 0.0 ;
	z= 0.0;

        angle = 0.0;
	tvPos = 0.0;
    }

};

#ifdef CROSS_SEV
class OccupiedIntersection : public CompactListItem<OccupiedIntersection>
{
    public:
    OccupiedIntersection(IntersectionItr is) : intersectionToCut(is) {};
    IntersectionItr intersectionToCut;
};
#endif

class TrajectoryVector
{
    public:

    Trajectory* trajectory;
    double pos;
    Direction dir,Cdir;

    TrajectoryVector(double p, Direction d, Trajectory* t) : pos(p), dir(d), trajectory(t) {}
    TrajectoryVector() {}

    TrajectoryVector flipped()
    {
        return TrajectoryVector(pos, dir == ::forward ? backward : ::forward  ,trajectory);
    }
};

struct OrderParameters
{
    double R,C;
    double Rdirector[3];
    vector<double> localOrder;
    vector<Vector3d> Sv;
};

class OrderParametersRaw
{
    public:

    double localL;
    double localLOpt;
    double Qxx;
    double Qxy;
    double Qxz;
    double Qyy;
    double Qyz;
    double Qzz;
    double localOrder;
    Vector3d Sv;

    OrderParametersRaw(void)
    {
        localL = 0.;
        localLOpt = 0.;
        Qxx = 0;
        Qxy = 0;
        Qxz = 0;
        Qyy = 0;
        Qyz = 0;
        Qzz = 0;

	Sv << 0.0,0.0,0.0;

        return;
    }

    double extractR(double director[3],string geometry)
    {
        double matrix[3][3] = {{Qxx,Qxy,Qxz},{Qxy,Qyy,Qyz},{Qxz,Qyz,Qzz}};
        double evecMat[3][3];
        double selEigenVal,eVal[3];
        eigen_decomposition(matrix, evecMat, eVal);

	if(geometry == "2D-plane_1_0_0")
	{
		int maxPos(0);
		if(eVal[1] > eVal[0])
		    maxPos = 1;
		if(eVal[2] > eVal[maxPos])
		    maxPos = 2;
		selEigenVal = fabs(eVal[maxPos]);
		for(int i=0 ; i<3 ; i++)
		director[i] = evecMat[i][maxPos];
	}
	else
	{
		int minPos(0);
		if(eVal[1] < eVal[0])
		    minPos = 1;
		if(eVal[2] < eVal[minPos])
		    minPos = 2;
		selEigenVal = fabs(eVal[minPos]);
		for(int i=0 ; i<3 ; i++)
		director[i] = evecMat[i][minPos];
	}

        return(selEigenVal);
    }
};

class Intersection
//this object describes the intersection of one trajectory with another, from one side.
//since it is created twice for every intersection, it is by far the most memory-intensive object
//and should be kept as small as possible
{
    public:
    int occupancy;
    Trajectory* otherTrajectory;
    IntersectionItr mirror;
    #ifdef CROSS_SEV
    OccupiedIntersection* occupiedListPtr;
    #endif
};

class PointATedge
{
    public:
    double x;
    double y;
    double z;

    int nextElement;
    PointATedge(double x1,double y1,double z1,int n1){x = x1;y = y1;z = z1;nextElement = n1;}
};

class Trajectory : public DLBaseItem<Trajectory>
//the trajectory is the basic geometrical object. Tips and segments associate with a trajectory
//and move/lie alongside it. Trajectory intersections determine the collision points
{
    friend class Region;
    public:

    const SurfaceVector base;
    const double length;

    vector<PointATedge>     endPoint;
    TrajectoryVector thisTr;

    // trajectoryVector of the connecting trajectory at the zero end
    TrajectoryVector prevTr;

    // trajectoryVector of the connecting trajectory at the far end
    TrajectoryVector nextTr;

    // cosine of the 3D angle with the previous trajectory (for edge catastrophes)
    double prevTrCosAngle;

    // cosine of the 3D angle with the next trajectory (for edge catastrophes)
    double nextTrCosAngle;

    // regular catastrophe value at the trajectory zero end
    double prevTrpCat;

    // regular catastrophe value at the trajectory far end
    double nextTrpCat;

    multimap<double, Intersection, std::less<double>, boost::fast_pool_allocator<std::pair<double, Intersection> > > intersections;

    // sorted list of all intersections
    IntersectionItr wallEnd()
    {
        // returns an iterator to the intersection that stands for the far wall
        return intersections.end();
    }
    IntersectionItr wallBegin()
    {
        // returns an iterator to the intersection that stands for the zero wall
        return intersections.begin();
    }

    explicit Trajectory(SurfaceVector,vector<PointATedge>,double);
    ~Trajectory();

    bool integrityCheck();

    // returns the connected trajectory in a given direction, and creates it if necessary
    TrajectoryVector nextTrajectory(Direction);

    // checks whether the trajectory can safely be removed - and does it
    void conditionalRemove();

    // list of pointers to associated segments
    list<Segment*> segments;

    // inserts a segment into the list
    TrjSegmentTag insertSegment(Segment*);

    // removes a segment from the list
    void removeSegment(TrjSegmentTag);

    // list of pointers to associated tips
    list<MTTip*> notificationList;

    // insert a tip for  notifications
    TrjMTTipTag registerForNotifications(MTTip*);

    // removes a tip from the list
    void unregisterForNotifications(TrjMTTipTag);

    // is called whenever an intersection is invalidated by trajectory removal
    void invalidateIntersection(IntersectionItr&);

    // is called when a new intersection is created
    void newIntersection(IntersectionItr&);

    // returns the difference sign of pos1-pos2, or itr1-itr2 if the first, cannot be determined accurately
    int differenceSign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2);

    // total length of trajectory (optical) that is covered by segments
    double coveredLength();

    // total length of segments on the trajectory
    double segmentLength();

    private:
    // avoid accidental (expensive) copying of Trajectory objects, by declaring private copy constructors without definitions
    Trajectory(const Trajectory&);
    Trajectory& operator=(const Trajectory&);
};


struct bendingOperator
{
    public:
        double cosTheta;
        double sinTheta;

        Vector3d axis;

        bendingOperator()
        {
            cosTheta = 0.0;
            sinTheta = 0.0;

            axis << 0.0,0.0,0.0;
        }
};

class Edge
{
public:
    int excludePoint;

    double pCat;
    double edgAngNorm;
    double xyRotationAngle;
    double edgAngle;
    double edgBendAngle;

    Vector2d midPoint;
    Vector2d b;
    Matrix2d A;
    vector<int> dir;
    Vector2d dir2D;
    vector<int> orientation;

    map<int,int> orientationMap;
    map<int,int> excludePointMap;

    Edge()
    {
        excludePoint = 0;

        pCat = 0.0;
        edgAngNorm = 0.0;
        xyRotationAngle = 0.0;
        edgAngle = 0.0;
	edgBendAngle = 0.0;

        midPoint << 0.0,0.0;

        A << 0.0,0.0,
             0.0,0.0;

        b << 0.0,0.0;

        orientation.push_back(0);
        orientation.push_back(0);
    }
};

class Region
//virtual base class for a geometry patch with a 2D coordinate system. Within a region,
//region functions are only called by geometry or trajectory member functions
{
    friend class Geometry;
    friend class Trajectory;

    public:

    // pointer to geometry
    Geometry* const geometry;

    //const RegionType type;
    double area;
    int regionId;
    double zOffset;
    double rotAngle;
    double periMeter;
    Vector2d midPoint;
    vector<Edge> side;
    map<int,int> sideMap;
    map<int,int> sideRevMap;
    Quaternion<double> Q;
    vector<Vector2d> Vertics;
    vector<Vector3d>orientation;
    int faceTag;
    int polyIntersectMark;
    vector<int> intersectEdg;

    // list of trajectories on the region
    DLList<Trajectory, TRAJECTORY_GRANULARITY> trajectories;

    // tip management functions
    list<MTTip*> growingPlusTipList;
    list<MTTip*> shrinkingPlusTipList;
    list<MTTip*> minusTipList;
    RegionMTTipTag registerOnRegion(MTTip*,TipType,MTType);
    void unregisterFromRegion(RegionMTTipTag,TipType,MTType);
    double totalLength;

    // time tag at which the length was last updated
    int previousUpdateTag;

    // updates the length to the current system time tag
    void updateRegionLength(bool forceUpdate = false);

    Region(Geometry* g, double a)
        :geometry(g), area(a),totalLength(0),
         previousUpdateTag(0),regionId(0),zOffset(0.0),rotAngle(0.0),periMeter(0.0)
         {
            midPoint << 0.0,0.0;
            for(int i=0;i<3;i++)
            {
                Vertics.push_back(Vector2d(0,0));
                side.push_back(Edge());
            }
        }
    virtual ~Region()
    {
        return;
    };
    double opticalLength();

    // vector computation functions (to be implemented in specializations)
    virtual void translateVector(SurfaceVector&, const double) =0;
    virtual void getTrajectoryCoordinates(SurfaceVector&, double&, vector<PointATedge>&,TrajectoryVector&) =0;
    virtual double intersectionAngle(Trajectory*, Trajectory*) =0;
    virtual void outputSnapshot(ostream&) = 0;
    virtual void outputOrderHeatMap(ostream&,vector<double>&,vector<Vector3d>&) = 0;

    protected:
    // trajectory management functions
    TrajectoryVector insertTrajectory(const SurfaceVector&);
    void removeTrajectory(Trajectory*);
    virtual void makeIntersectionList(Trajectory*) =0;
    virtual SurfaceVector randomSurfaceVector() =0;

    private:
    // avoid accidental (expensive) copying of Region objects, by declaring private copy constructors without definitions
    Region(const Region&);
    Region& operator=(const Region&);
};

class elementList
{
    public:
    int midPoint;
    double length;
    vector<int> nodes;
    vector<int> sharedElement;
    vector<int> sharedEdge;
    double edge3DAngle;
    Vector3d intersectionByPlane;

    elementList()
    {
        midPoint = 0;
        length = 0.0;
	edge3DAngle = 0.0;
        intersectionByPlane << 0.0,0.0,0.0;

        for(int i=0;i<2;i++)
        {
            nodes.push_back(0);
            sharedElement.push_back(0);
            sharedEdge.push_back(0);
        }
    }
};

class Geometry
//virtual base class for the various types of geometries
{
    friend class System;

    public:
    //const GeometryType type;
    double area;
    double areainMesh;
    System* const system;

    vector<elementList> ppbEdgeList;
    vector<Vertics> globalVertex;


    double areaPPB;
    vector<Vertics> ppb;

    double patchArea[3];
    vector<int> RegionsIndex[3];

    int elementMax;
    Vector3d objectCM,objectPA,nucleousPosition;
    vector<Region*> regions;
    Geometry(System* s,double a) : area(a),system(s) {for(int i=0;i<3;i++)patchArea[i]=0.0;};

    // virtual destructor implies virtual class
    virtual ~Geometry()
    {
        return;
    };
    bool integrityCheck();
    double opticalLength();
    int trajectoryCount();
    void callTranslator(SurfaceVector&,int,int);

    virtual void getOrderParameters(OrderParameters&) = 0;
    virtual void outputSnapshot(ostream&) = 0;
    virtual void outputOrderHeatMap(ostream&,vector<double>&,vector<Vector3d>&) = 0;

    SurfaceVector randomSurfaceVector();
    TrajectoryVector createTrajectory(const SurfaceVector&);

    // creates a trajectory at a given surface vector and returns a trajectory vector
    TrajectoryVector createAndLinkTrajectory(const SurfaceVector&, Trajectory*, Direction, double,double);

    // only called by Trajectory (friend class)
    private:
    friend TrajectoryVector Trajectory::nextTrajectory(Direction);
    virtual TrajectoryVector extendTrajectory(Trajectory*, Direction) =0;
};

class DeterministicQueue
//'Smart' queue object that contains logic to transform to and from the system time.
{
    //public:
    friend class System;
    priority_queue<DeterministicEvent> queue;	// the queue itself
    double currentBase;							// the current 'distance'
    double valueCache[POSITION_CACHE_SIZE];		// cache of previous 'distance' values
    double (System::* const distanceTimeConversionFunction)(double);
    // pointer to a function that converts a distance to a time
    double (System::* const timeDistanceConversionFunction)(double);
    // pointer to a function that converts a time to a distance

    public:
    System* const 	system;						// pointer to containing system
    DeterministicQueue(System*, double (System::*dtFunc)(double), double (System::*tdFunc)(double));
    double currentPos()
    {
        return currentBase;    // returns the current 'distance'
    }
    void advanceTime(double);					// advances the parameters to the current system time
    void storeTime(int);						// store 'currentBase' in the cache at a given tag position
    double progression(int cachePos)
    {
        return currentBase - valueCache[cachePos];
    }
    // returns the progress of 'currentBase' relative to a given time tag
    double firstEventTime();					// returns the system time of the first scheduled event
    DeterministicEvent pop();					// returns and removes the first scheduled event
    bool empty()
    {
        return queue.empty();    // checks whether the queue is empty
    }
    void flush();								// empties the queue and resets currentBase
    void pushGlobal(double, GlobalEventType);	// pushes a global event onto the queue at a certain distance
    EventTrackingTag pushDeterministic(double, EventDescriptorIndex);
    // pushes a deterministic (MT) event onto the queue at a certain distance
};

class EventDescriptor
{
    //The EventDescriptor object is the microtubule-side interface to the event queues.
    //It stores information about the upcoming event, and is used as an entry-point for the
    //event handler when this event is executed.
    //Every microtubule has three of these objects: one on every tip and one to keep track of
    //'ev_disappear' events, when the two tips annihilate each other.

    public:

    // index of this object in the EventDescriptorMap (necessary for removal)
    EventDescriptorIndex index;
    Microtubule* const 		mt;

    // type of queued event
    DeterministicEventType 	type;

    // associated tracking tag (for validity testing)
    EventTrackingTag 		tag;

    // inverse velocity (with respect to the event clock defined in the queue)
    double distanceScaleFactor;
    DeterministicQueue* 	queue;
    EventDescriptor(Microtubule*, DeterministicQueue*, double);
    ~EventDescriptor();

    // resets the queue and velocity of the event timer
    void reinitialize(DeterministicQueue*, double) ;

    // pushes a deterministic event at a certain distance onto the queue
    void pushOnQueue(double, DeterministicEventType);

    // invalidates the queued event
    void clear();
};

class DeterministicEvent
{
    public:
    double eventTimeDist;

    // -1 indicates a global event
    EventDescriptorIndex infoIdx;

    union
    {
        EventTrackingTag tag;
        GlobalEventType	global_type;
    };

    // define the < operator for automatic sorting of events in the event queue (nearest first)
    bool operator<(const DeterministicEvent& ev2) const
    {
        return eventTimeDist > ev2.eventTimeDist;
    }
};

class MTTip
{
    public:
    Microtubule* mt;
    Trajectory* trajectory;

    TrjMTTipTag notificationTag;
    RegionMTTipTag regionTag;
    double velocity;
    Direction dir;
    IntersectionItr	nextCollision;
    double nextEventPos;
    EventDescriptor event;
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
    // avoid accidental (expensive) copying of Trajectory objects, by declaring private copy constructors without definitions
    MTTip(const MTTip&);
    MTTip& operator=(const MTTip&);

};

class Segment : public DLBaseItem<Segment>
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

    // constructs a segment as part of a microtubule at a specified vector location
    Segment(Microtubule*, TrajectoryVector&);

    ~Segment();
    double length()
    {
        return abs(end-start);
    }
    bool isLastInMT();
    bool isFirstInMT();
    bool crossesIntersection(IntersectionItr& is) ;
};

class Microtubule : public DLBaseItem<Microtubule>
{
    public:
    System* const system;
    double nucleationTime;

    // time tag at which the positions and lengths were last updated
    int previousUpdateTag;

    // type of microtubule (growing, shrinking)
    MTType type;
    MTTip plus;
    MTTip minus;
    DLList<Segment>	segments;

    // event descriptor for events that are not directly associated to a single tip (only disappearance)
    EventDescriptor disappearEvent;

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
    void sever(Segment*,double);
    void severAtCross(IntersectionItr is, Segment* cutSeg);
    void splitSegmentAtTrajPos(double cutPos, Segment* cutSeg);

   // translates random position at MT to corresponding position at Segment
    void translatePositionMT2Segment(double& cutPos, Segment*& cutSeg);

    private:
    // avoid accidental (disastrous) copying of Trajectory objects, by declaring private copy constructors without definitions
    Microtubule(const Microtubule&);
    Microtubule& operator=(const Microtubule&);
};

class Measurement
{
    public:
    double time;
    double lengthDensity;
    double opticalDensity;
    double averageLength;
    int numberOfMTs;
    double segmentsPerMT;
    OrderParameters order;
    int growingNumber;
    int shrinkingNumber;
    int segments;
    int trajectories;
    int zipperCount;
    int crossoverCount;
    int inducedCatastropheCount;
    int validDEventCount;
    int invalidDEventCount;
    int sEventCount;
    int lengthSeveringCount;
    int intersectionSeveringCount;
    int occupiedIntersectionCount;
    double G_effAdjust;
};

class Parameters
{
    public:
    System* system;
    double vPlus;
    double vMin;
    double vTM;
    double kSev;
    double kCross;
    double kCat;
    double kRes;
    double kNuc;
    double PPBformingTime;
    double poolDensity;
    int treadmillingEnabled;
    int severingEnabled;
    int crossSeveringEnabled;
    int crossSeveringTop;
    double crossSeveringStartAngle;
    int PPB;
    int restrictedPool;
    int edgeCatastropheEnabled;
    int edgeCatastropheSmooth;
    double pCatSpecialEdgeMax;
    double pCatRegularEdgeMax;
    int edgNumber;
    int faceNumber;
    map<string,double> edgCatMap;
    map<string,double> faceCatMap;
    vector<double> RegionKcatMultiplier;
    NucleationType nucleationType;
    int zipperingEnabled;
    int catastrophesEnabled;
    int proportionalCatastrophes;
    double inducedCatastropheFraction;
    double zipFraction;
    double PPBkNucFraction;
    double catStartAngle;
    double magicAngle;
    InteractionType interactionType;
    BundleType bundleType;
    string geomParam;
    unsigned long seed;
    double stopTime;
    double measurementInterval;
    double wallClockLimit;
    double memoryLimit;
    string inputDir;
    string outputDir;
    string movieDir;
    int createSubdir;
    double newParameterReadInterval;
    string newParameterFile;
    int movieEnabled;
    int showOutput;
    int showMesh;
    string geometry;
    double movieFrameInterval;
    double geometryScaleFactor;

    double c0calc;
    double x0calc;
    double z0calc;

    public:
    Parameters(System*);
    void initialize(const char *);
    bool reinitialize(const char *);
    bool readFromFile(const char*, bool);
    bool writeToFile();
    void verifyParameters();
    bool calcTheoryParameters();
};

class System
//Master object defining an interacting mt system
{
    public:
    Geometry* geometry;
    DLList<Microtubule,MICROTUBULE_GRANULARITY> growing_mts;
    DLList<Microtubule,MICROTUBULE_GRANULARITY> shrinking_mts;

    #ifdef CROSS_SEV
    CompactList<OccupiedIntersection,OCCUPIED_INTERSECTION_GRANULARITY> OccupiedIntersectionList;
    #endif

    vector<int> growingTipsReg;

    Parameters p;

    // random generator
    MTRand randomGen;

    // current system time with respect to offset
    double systemTime;

    // current system time offset
    double systemTimeOffset;

    int currentTimeTag;
    double totalLength;
    int memUsage;
    bool stopSignal;

    time_t wallClockStartTime;
    int totalSEventCount;
    int totalValidDEventCount;
    int totalInvalidDEventCount;
    int totalLengthSeveringCount;
    int totalIntersectionSeveringCount;
    int boundaryCrossingCount;
    int totalZipperCount;
    int totalCrossoverCount;
    int totalInducedCatastropheCount;

    int countSegments;
    int countTrajectories;
    int countIntersections;
    int estimateMemoryFootprint();

    double nextStochasticEventTime;
    StochasticEventType nextStochasticEventType;

    // queue object for time-defined events
    DeterministicQueue timeQueue;

    // queue object for vPlus-defined events (regulated by tubulin pool size)
    DeterministicQueue vPlusQueue;

    double identity(double i)
    {
        return i;
    }

    // distance to time conversion for growing tips
    double vPlusToTime(double);

    // time to distance conversion for growing tips
    double timeToVPlus(double);

    double nextStatusEventTime;
    double nextSnapshotEventTime;
    double nextParameterEventTime;
    double nextPPBEventTime;

    // current tag for event descriptor
    EventDescriptorIndex EventDescriptorID;

    // structure that maps EventDescriptorIDs onto EventDescriptor objects
    map<EventDescriptorIndex, EventDescriptor*> EventDescriptorMap;

    // returns a new event descriptor tag
    EventDescriptorIndex registerEventDescriptor(EventDescriptor*);

    // invalidates an event descriptor tag
    void unregisterEventDescriptor(EventDescriptorIndex);

    // returns a pointer to an EventDescriptor through its tag
    EventDescriptor* getEventDescriptor(EventDescriptorIndex);

    // the current event tracking tag
    EventTrackingTag eventID;

    // returns a new event tracking tag
    EventTrackingTag getEventTag()
    {
        return eventID++;
    }

    System(char*);
    ~System();
    bool integrityCheck();
    void flushAndReload(bool = true);
    void advanceTime(double nextTime);

    // updates all positions in the system to the system time
    void updateAll(bool forceUpdate = false);
    void emergencyBreak();

    void run(double,string&,string&);
    void nextEvent(void);
    void handleGlobalEvent(DeterministicEvent&);
    void determineStochasticEvent();
    void handleNucleationEvent();
    void handleSeveringEvent();

    #ifdef CROSS_SEV
    void handleSeveringAtCrossEvent();
    #endif

    void handleRescueEvent();
    void handleCatastropheEvent();
    void randomPositionOnMicrotubule(double& cutLength, Segment*& cutSeg);

    // returns the result of a collision at a specified angle, and bundle occupancies
    CollisionType collisionResult(double, int, int);

    // returns the probabilities for induced catastrophes and zippering, as a function of the angle
    void collisionProbabilities(double, double&, double&);


    // binomial table. Necessary for bundle collision type: multi-collisions
    double binomialTable[MAXBINOM][MAXBINOM];
    void makeBinomialTable(void);

    double multiPcross(int Npar, int Ncoll, double xSingle, double zSingle);
    double multiPzip(int Npar, int Ncoll, double xSingle, double zSingle);

    ofstream parameterFile;
    ofstream movieFile;
    ofstream heatMapFile;
    ofstream measurementFile;
    ofstream orderDirectorFile;
    deque<Measurement> measurementHistory;
    bool dataSaved;
    void initializeOutput();
    void performMeasurement();
    void writeMeasurementsToFile(int = 0);
    void closeFiles();

    #ifdef CROSS_SEV
    // removes pointers from is and its mirror
    void removeOccupiedIntersection(Intersection& is);

    // creates pointers for is and its mirror [occupiedIntersection contains only one pointer: to the one "on top" (to be cut)]
    void addOccupiedIntersection(IntersectionItr is);
    #endif

    private:
    // avoid accidental (expensive) copying of Trajectory objects, by declaring private copy constructors without definitions
    System(const System&);
    System& operator=(const System&);
};

ostream& operator<<(ostream&, const Measurement);
void writeMeasurementDescriptors(ostream&);

class TriMeshGeometry: public Geometry
{
    public:
    TriMeshGeometry(System*);
    TrajectoryVector extendTrajectory(Trajectory*, Direction);
    void getOrderParameters(OrderParameters&);
    void outputSnapshot(ostream&);
    void outputOrderHeatMap(ostream&,vector<double>&,vector<Vector3d>&);
};

class Cartesian : public Region
{
    public:
    Cartesian(Geometry*, double);

    virtual SurfaceVector randomSurfaceVector() = 0;
    virtual void getTrajectoryCoordinates(SurfaceVector&, double&,vector<PointATedge>&,TrajectoryVector&) = 0;
    void makeIntersectionList(Trajectory*);
    void translateVector(SurfaceVector&, const double);
    double intersectionAngle(Trajectory*, Trajectory*);

    void getOrderParameters(OrderParameters&);
    void getOrderParametersRawFlat(OrderParametersRaw&,vector<Vector3d>&,double);
    void getOrderParametersRawCylinder(OrderParametersRaw&, double, double, double);

    void outputSnapshot(ostream&);
    void outputSnapshotOffset(ostream&,double,double);

    void outputOrderHeatMap(ostream&,vector<double>&,vector<Vector3d>&);
};

class Triangle : public Cartesian
{
    public:
    Triangle(double, Geometry*);

    SurfaceVector randomSurfaceVector();
    void getTrajectoryCoordinates(SurfaceVector&, double&,vector<PointATedge>&,TrajectoryVector&);

};

#ifndef NO_INLINE
#include "inline.cpp"
#endif

#endif
//CORTICALSIM_H_



