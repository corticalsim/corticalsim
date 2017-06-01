/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef CORTICALSIM_H_
#define CORTICALSIM_H_

/* corticalSim - Cortical microtubule simulator
 * Simon Tindemans and Eva Deinum, FOM Institute AMOLF, 2008-2017
 * 
 * IMPORTANT: pitfalls to watch out for when extending the code
 * * A right handed 3D coordinate system is used, based on screen coordinates. X increased to the right, Y increases
 *   upward and Z increases toward the viewer.
 * * A physical trajectory on a region has two equivalent orientations. However, we must assign an orientation for the
 *   trajectory coordinates and their transitions. The mapping between surface vectors and trajectory vectors should be 
 *   1 on 1, otherwise inconsistencies can occur. [corrupt trajectory links, etc]
 * * The microtubule is only assumed to grow in length at the plus end. This is not an MTTip property, but a 
 *   property of the Microtubule functions. Two-sided growth can be enabled, but needs plus/minus distinctions
 *   for the events. Also, tips with velocity zero (usually non-treadmilling minus ends) are
 *   assumed to be in a 'shrinking' state.
 * * The program assumes a stable sort order for the intersection lists. This is true in practice (in all common
 *   implementations as far as I could ascertain), but not mandated by the C++ standard.
 * * In microtubule.cpp and system.cpp, compiler warnings are (should be) issued for using the 'this' pointer within
 *   the initializer list, because the 'this' pointer cannot be used until the initialization has completed. The 
 *   current use is ok, because the pointer is only used to store the address of the parent object. NOTE: when making 
 *   changes to the code, take care NOT to use this pointer within the constructors of other objects!

 TO DO:
 - remove question marks and exclamations in intersection ordering code (made by eva?)
 */


#define PROGRAM_VERSION "1.26.1"

//#define NO_INLINE				// if defined, avoids inlining the functions in inline.cpp. This 
								// carries a performance penalty, but may aid debugging
//#define NO_INTEGRITY_CHECK	// if defined, disables the system integrity checks performed
								// at every measurement. This may speed up execution if you perform 
								// measurements at a high rate, but makes the code less aware of
								// potential bugs. Not recommended.

// various debug flags which can be enabled at compile time
//#define DBG_ASSERT
//#define DBG_ACID_TEST
//#define DBG_QMGMT
//#define DBG_CONSTRUCT
//#define DBG_GEOMETRY
//#define DBG_MTS
//#define DBG_EVENT
//#define DBG_SYSTEM
//#define DBG_INTERACTION_TEST
	
// disable specific compiler warnings
#define NOMINMAX					// to prevent problems with numeric_limits on VC++

#define BOOST_DISABLE_THREADS      // permit header-only compilation of Boost::pool
									// recent versions support threading which requires compilation

//#pragma warning(disable:981)
//#pragma warning(disable:522)
//#pragma warning(disable:383)

#include <vector>
#include <list>
#include <map>
#include <queue>
#include <boost/pool/pool_alloc.hpp>
#include <boost/lexical_cast.hpp>  //Wat doet dit?
#include <string>
#include <cstring>
#include <cmath>
#include <fstream>
#include "MersenneTwister.h"
#include "eig3.h"
#include "DLList.h"
#include "CompactList.h"


using namespace std;

// numeric parameters
const double PI = 3.141592653589793;
const double ZERO_CUTOFF = 100000*numeric_limits<double>::epsilon();   // approx 10E-10;  
const double VERY_LARGE = 10E100;

const int MAXBINOM = 100;

// memory management
const int MICROTUBULE_GRANULARITY = 256;
const int SEGMENT_GRANULARITY = 1024;
const int TRAJECTORY_GRANULARITY = 512;
const int INTERSECTION_GRANULARITY = 4096;
const int EVENT_GRANULARITY = 1024;
const int OCCUPIED_INTERSECTION_GRANULARITY = 256;

const int MAX_HISTORY_SIZE = 256; // maximum size of measurement cache (increase to minimize disk access)
const int MIN_HISTORY_SIZE = 64;  // minimum size of measurement cache (minimum may be needed for trending) 

// simulation parameters
const int POSITION_CACHE_SIZE = 32768;
const int CLOCK_POLLING_INTERVAL = 10000;
const int MEMORY_POLLING_INTERVAL = 200;
const int QUEUE_FLUSH_INTERVAL = 1000000;			//note that this is related to the total number of events...

const int NUCLEATION_DISCRETIZATION_STEPS = 512;



typedef enum {g_periodic, g_grid, g_wormhole, g_box, g_cylinder, g_gridcylinder, g_spherocylinder, g_sphere, g_pancake, g_COUNT_LAST} GeometryType;
extern string GeometryTypeText[]; // = {"periodic", "grid", "wormhole", "box", "cylinder", "gridcylinder", "spherocylinder", "sphere", "pancake"};

typedef enum {r_rectangle, r_disc, r_dome} RegionType;
extern string RegionTypeText[]; // = {"rectangle", "disc", "dome"};

typedef enum {forward=1, backward=-1} Direction;
typedef enum {s_top, s_left, s_bottom, s_right} RectangleSide;

typedef enum {ev_none, ev_wall=1, ev_collision, ev_backtrack, ev_end_of_segment, ev_disappear} DeterministicEventType;
typedef enum {status=1, snapshot, stop, parameter_change} GlobalEventType;
typedef enum {catastrophe=1, rescue, katanin, severingAtCross, nucleation, preSeededNucleation} StochasticEventType;
typedef enum {ct_zipper=1, ct_crossover, ct_inducedCatastrophe} CollisionType;

//typedef enum {nuc_isotropic, nuc_biased, nuc_discreteAngles, nuc_chanLloyd, nuc_chanLloydRandomPosition, nuc_chanLloydIsotropic, nuc_ellipse,  nuc_COUNT_LAST} NucleationType;	
//extern string NucleationTypeText[]; // = {"isotropic", "biased", "discreteAngles", "chanLloyd", "chanLloydRandomPosition", "chanLloydIsotropic", "ellipse" };
typedef enum {nuc_isotropic, nuc_biased, nuc_discreteAngles, nuc_ellipse,  nuc_COUNT_LAST} NucleationType;	
extern string NucleationTypeText[]; // = {"isotropic", "biased", "discreteAngles", "ellipse" };

typedef enum {int_zipFirst,  int_catFirst, int_minimalFourier, int_COUNT_LAST} InteractionType;
extern string InteractionTypeText[];// = {"zipFirst", "catFirst", "minimalFourier" };

typedef enum {bdl_simple,  bdl_sticky, bdl_noZip, bdl_multiCollision, bdl_Ncollision, bdl_COUNT_LAST} BundleType;
extern string BundleTypeText[];// = {"simple", "sticky", "noZip", "multiCollision", "Ncollision" };

typedef enum {t_minus, t_plus} TipType;
typedef enum {mt_growing, mt_shrinking} MTType;

typedef enum {graphics_mt, graphics_trajectory, graphics_plus, graphics_minus} GraphicsType;

typedef enum { r_accounting_normal, r_accounting_dontcount } RegionAccountingType;
typedef enum { r_param_normal, r_param_modified } RegionParameterType;

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
class OccupiedIntersection;

typedef multimap<double, Intersection, std::less<double>, boost::fast_pool_allocator<std::pair<double, Intersection> > >::iterator IntersectionItr;
typedef int EventDescriptorIndex;	// TO DO: turn into unsigned, and replace -1 by ::max()
typedef int EventTrackingTag;		// TO DO: turn into unsigned, and replace -1 by ::max()

typedef list<Segment*>::iterator TrjSegmentTag;
typedef list<MTTip*>::iterator TrjMTTipTag;
typedef list<MTTip*>::iterator RegionMTTipTag;

class Coord2D {
public:
	union { 
	double x;
	double a;
	};
	union {
	double y;
	double b;
	};
	Coord2D(double v1, double v2) : x(v1), y(v2) {};
};

struct SurfaceVector {
	/* defines a tangent vector on the surface of a geometry.
	 * 
	 * Note that the names x, y and angle are chosen to reflect their meaning for 
	 * for 2D cartesian surfaces, but their meaning may differ depending on the region type
	 */
	double x;
	double y;
	double angle;
	Region* region;
};

class OccupiedIntersection : public CompactListItem<OccupiedIntersection>
{
public:
	OccupiedIntersection(IntersectionItr is) : intersectionToCut(is) {}; 
	IntersectionItr intersectionToCut;
};


class TrajectoryVector {
public:
	Trajectory* trajectory;		// pointer to trajectory
	double pos;					// position along the trajectory
	Direction dir;				// direction, relative to the trajectory

	// two different constructors to allow for flexibility in object creation
	TrajectoryVector(double p, Direction d, Trajectory* t) : pos(p), dir(d), trajectory(t) {}
	TrajectoryVector() {}
	
	TrajectoryVector flipped() { return TrajectoryVector(pos, dir == ::forward ? backward : ::forward  ,trajectory); }
								// returns a trajectory vector with the opposite direction
};



/*************************** Analysis tools ********************************/


struct OrderParameters
{
	double S2;
	double S2angle;
	double S4;
	double S4angle;
	double S2Opt;
	double S2angleOpt;
	double S4Opt;
	double S4angleOpt;

	double R;
	double Rdirector[3];
};

class OrderParametersRaw
{
public:
	double si2;
	double si4;
	double co2;
	double co4;
	double localL;
	double si2Opt;
	double si4Opt;
	double co2Opt;
	double co4Opt;
	double localLOpt;

	double Qxx;
	double Qxy;
	double Qxz;
	double Qyy;
	double Qyz;
	double Qzz;

	double isoWeights[3];

	OrderParametersRaw(void)
	{
		si2 = 0.;
		si4 = 0.;
		co2 = 0.;
		co4 = 0.;
		localL = 0.;
		si2Opt = 0.;
		si4Opt = 0.;
		co2Opt = 0.;
		co4Opt = 0.;
		localLOpt = 0.;
		Qxx = 0;
		Qxy = 0;
		Qxz = 0;
		Qyy = 0;
		Qyz = 0;
		Qzz = 0;

		isoWeights[0] = 0;
		isoWeights[1] = 0;
		isoWeights[2] = 0;
		return;
	}

	double extractR(double director[3])
	{
		double matrix[3][3] = {{Qxx-isoWeights[0],Qxy,Qxz},{Qxy,Qyy-isoWeights[1],Qyz},{Qxz,Qyz,Qzz-isoWeights[2]}};

		// assume sum(isoWeights)=1
/*		for (int i=0 ; i<3 ; i++)
			for (int j=0 ; j<3 ; j++)
				matrix[i][j] /= sqrt(1.0 - isoWeights[i])*sqrt(1.0 - isoWeights[j]);
*/
		double evecMat[3][3];
		double eVal[3];
		int minPos;

		eigen_decomposition(matrix, evecMat, eVal);
		minPos = 0;
		if (eVal[1] < eVal[0])
			minPos = 1;
		if (eVal[2] < eVal[minPos])
			minPos = 2;
		for (int i=0 ; i<3 ; i++)
			director[i] = evecMat[i][minPos];
		return (-eVal[minPos]/max(isoWeights[0],max(isoWeights[1],isoWeights[2])));
	}

};


/************************** Geometry/regions/trajectories/intersections ***************************/


class Intersection
/*
 * The intersection object describes the intersection of one trajectory with another, from one side.
 * Since it is created twice for every intersection, it is by far the most memory-intensive object
 * and should be kept as small as possible.
 */
{
public:
	int occupancy;					// number of MTs that are crossing the intersection along this trajectory
	Trajectory* otherTrajectory;	// pointer to the other trajectory
	IntersectionItr mirror;			// pointer to the mirroring Intersection structure of the other trajectory
	// TODO: maybe add corner count
	OccupiedIntersection* occupiedListPtr;
};

class Trajectory : public DLBaseItem<Trajectory>
/*
 * The trajectory is the basic geometrical object. Tips and segments associate with a trajectory
 * and move/lie alongside it. Trajectory intersections determine the collision points.
 * 
 */
{
	friend class Region;

public:
	const SurfaceVector		base;			// coordinates and angle of the base of the trajectory
	const double 			length;			// total length of the trajectory, from one boundary to another
	TrajectoryVector		prevTr;			// continuation of the Trajectory at the base
	TrajectoryVector		nextTr;			// continuation of the Trajectory at the far end

	double					prevTrCosAngle;	// Cosine of the 3D angle with the previous trajectory (for edge catastrophes)
	double					nextTrCosAngle; // Cosine of the 3D angle with the next trajectory (for edge catastrophes)

	multimap<double, Intersection, std::less<double>, boost::fast_pool_allocator<std::pair<double, Intersection> > > intersections;
																	// sorted list of all intersections
	IntersectionItr wallEnd() { return intersections.end(); }		// returns an iterator to the intersection that stands for the far wall
	IntersectionItr wallBegin() { return intersections.begin(); }	// returns an iterator to the intersection that stands for the zero wall

	explicit Trajectory(SurfaceVector, double);	
	~Trajectory();
	bool integrityCheck();
	TrajectoryVector nextTrajectory(Direction);						// returns the connected trajectory in a given direction, and creates it if necessary
	void conditionalRemove();										// checks whether the trajectory can safely be removed - and does it

	// segment management functions
	list<Segment*> segments;										// list of pointers to associated segments
	TrjSegmentTag insertSegment(Segment*);							// inserts a segment into the list
	void removeSegment(TrjSegmentTag);								// removes a segment from the list

	// tip management functions
	list<MTTip*> notificationList;									// list of pointers to associated tips
	TrjMTTipTag registerForNotifications(MTTip*);					// registers a tip for insertion/deletion notifications
	void unregisterForNotifications(TrjMTTipTag);					// removes a tip from the list

	// taking care of insertion/removal of new intersections
	void invalidateIntersection(IntersectionItr&);					// is called whenever an intersection is invalidated by trajectory removal
	void newIntersection(IntersectionItr&);							// is called when a new intersection is created
	int differenceSign(IntersectionItr itr1, double pos1, IntersectionItr itr2, double pos2);
																	// returns the difference sign of pos1-pos2, or itr1-itr2 if the first
																	// cannot be determined accurately.

	// analysis functions
	double coveredLength();											// total length of trajectory that is covered by segments
	double segmentLength();											// total length of segments on the trajectory

private:
	// avoid accidental (expensive) copying of Trajectory objects, by declaring private copy constructors without definitions
	Trajectory(const Trajectory&);
	Trajectory& operator=(const Trajectory&);
};









class Region
/* 
 * Virtual base class for a geometry patch with a 2D coordinate system. Within a region,
 * trajectories may only intersect once.
 * 
 * Region functions are only called by geometry or trajectory member functions.
 */
{
	friend class Geometry;
	friend class Trajectory;

public:
	Geometry* const 	geometry;			// pointer to containing geometry
	const int geometryRegionIndex;			// region index corresponding to this region within the geometry.
	const RegionType	type;				// type of region (rectangle, disc, hemisphere)
	const double 		area;				// area of the region
	const RegionAccountingType accountingType;	// whether segments in this region are part of the 'main' region or not
	const RegionParameterType parameterType;

	DLList<Trajectory, TRAJECTORY_GRANULARITY> trajectories;
											// list of trajectories on the region

	// tip management functions
	list<MTTip*> growingPlusTipList;									// list of pointers to growing tips
	list<MTTip*> shrinkingPlusTipList;									// list of pointers to shrinking tips
	list<MTTip*> minusTipList;									// list of pointers to treadmilling tips
	RegionMTTipTag registerOnRegion(MTTip*,TipType,MTType);					// registers a tip for insertion/deletion notifications
	void unregisterFromRegion(RegionMTTipTag,TipType,MTType);					// removes a tip from the list
	double totalLength;

	int previousUpdateTag;			// time tag at which the length was last updated
	void updateRegionLength(bool forceUpdate = false);			// updates the length to the current system time tag



	Region(RegionType t, Geometry* g, int index, double a, RegionAccountingType acc, RegionParameterType par) 
		: type(t), geometry(g), geometryRegionIndex(index), area(a), accountingType(acc), parameterType(par), totalLength(0), previousUpdateTag(0) {};
	virtual ~Region() {return;};


	double opticalLength();

	// vector computation functions (to be implemented in specializations)
	virtual void translateVector(SurfaceVector&, const double) =0;
	virtual void getTrajectoryCoordinates(SurfaceVector&, double&, TrajectoryVector&) =0;
		// NOTE: it is crucial that implementations return the same (base) SurfaceVector for any SurfaceVector within the trajectory
	virtual double intersectionAngle(Trajectory*, Trajectory*) =0;
	virtual double calculateHistogram(double[], bool optical = false, bool=false) = 0;
	virtual void outputSnapshot(ostream&) = 0;

protected:
	// trajectory management functions
	TrajectoryVector insertTrajectory(const SurfaceVector&, bool scanExisting = false);	// inserts a trajectory at a given surfact vector
	void removeTrajectory(Trajectory*);							// removes a trajectory from the region
	virtual void makeIntersectionList(Trajectory*) =0;
	virtual SurfaceVector randomSurfaceVector() =0;				// creates a random surface vector

private:
	// avoid accidental (expensive) copying of Region objects, by declaring private copy constructors without definitions
	Region(const Region&);
	Region& operator=(const Region&);
};











class Geometry
/* 
 * Virtual base class for the various types of geometries
 */
{
	friend class System;

public:	
	const GeometryType 	type;					// the geometry type (periodic, sphere, etc.)
	const double 		area;					// total area of all surfaces in the geometry
	System* const		system;					// pointer to the containing system
//protected:
	vector<Region*> regions;					// list of the Region-derived objects that make up the geometry

public:
	Geometry(GeometryType g, System* s, double a) : type(g), area(a), system(s) {};
	// virtual destructor implies virtual class
	virtual ~Geometry(){ return; };
	bool integrityCheck();						// performs an integrity check of the trajectory network
	double opticalLength();
	int trajectoryCount();
	virtual void getOrderParameters(OrderParameters&) = 0;
	virtual double calculateHistogram(double[],bool=false) = 0;
	virtual void outputSnapshot(ostream&) = 0;

public:
	SurfaceVector randomSurfaceVector();		// returns a random surface vector
	TrajectoryVector createTrajectory(const SurfaceVector&, bool scanExisting = false);
												// creates a trajectory at a given surface vector and returns a trajectory vector
	TrajectoryVector createAndLinkTrajectory(const SurfaceVector&, Trajectory*, Direction, double);

private:	// only called by Trajectory (friend)
	friend TrajectoryVector Trajectory::nextTrajectory(Direction);
	virtual TrajectoryVector extendTrajectory(Trajectory*, Direction) =0;
};



/********************* control/event-related functions ****************************/

class DeterministicQueue
/*
 * 'Smart' queue object that contains logic to transform to and from the system time.
 * 
 */
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

	double currentPos() { return currentBase; }	// returns the current 'distance'
	void advanceTime(double);					// advances the parameters to the current system time
	void storeTime(int);						// store 'currentBase' in the cache at a given tag position
	double progression(int cachePos) { return currentBase - valueCache[cachePos]; }
												// returns the progress of 'currentBase' relative to a given time tag
	double firstEventTime();					// returns the system time of the first scheduled event
	DeterministicEvent pop();					// returns and removes the first scheduled event
	bool empty() { return queue.empty(); }		// checks whether the queue is empty
	void flush();								// empties the queue and resets currentBase

	void pushGlobal(double, GlobalEventType);	// pushes a global event onto the queue at a certain distance
	EventTrackingTag pushDeterministic(double, EventDescriptorIndex);
												// pushes a deterministic (MT) event onto the queue at a certain distance
};








class EventDescriptor
{
/* The EventDescriptor object is the microtubule-side interface to the event queues. 
 * It stores information about the upcoming event, and is used as an entry-point for the
 * event handler when this event is executed. 
 * 
 * Every microtubule has three of these objects: one on every tip and one to keep track of 
 * 'ev_disappear' events, when the two tips annihilate each other.
 */
public:
	EventDescriptorIndex index;		// index of this object in the EventDescriptorMap (necessary for removal)
public:
	Microtubule* const 		mt;		// pointer to the parent MT
	DeterministicEventType 	type;	// type of queued event
	EventTrackingTag 		tag;	// associated tracking tag (for validity testing)
	double distanceScaleFactor;		// inverse velocity (with respect to the event clock defined in the queue)
	DeterministicQueue* 	queue;	// pointer to the queue that is used

public:
	EventDescriptor(Microtubule*, DeterministicQueue*, double);
	~EventDescriptor();
	void reinitialize(DeterministicQueue*, double) ;		// resets the queue and velocity of the event timer

	void pushOnQueue(double, DeterministicEventType);		// pushes a deterministic event at a certain distance onto the queue
	void clear();											// invalidates the queued event
};









class DeterministicEvent {
public:
	double eventTimeDist;
	EventDescriptorIndex infoIdx;			// -1 indicates a global event
	union {
		EventTrackingTag tag;		// maybe enter tracking number automatically in constructor
		GlobalEventType	global_type;
	};

	// define the < operator for automatic sorting of events in the event queue (nearest first)
	bool operator<(const DeterministicEvent& ev2) const { return eventTimeDist > ev2.eventTimeDist; }
};



/************************** Bio stuff: MTs, segments, tips ******************************/






class MTTip
{
public:
	Microtubule* mt;				// pointer to the containing microtubule
	Trajectory* trajectory;			// pointer to the containing trajectory
	TrjMTTipTag notificationTag;	// Token for the tip-trajectory link
	RegionMTTipTag regionTag;		// Token for the tip-region link

	double velocity;				// velocity of the tip, with respect to the microtubule (positive is polymerizing) and the event clock
	Direction dir;					// polymerization direction with respect to trajectory (opposite for plus and minus tips on a single segment!)

//private:
	IntersectionItr	nextCollision;	// iterator to the next intersection
	double nextEventPos; 			// position of the next traversal/collision event on the trajectory
	EventDescriptor event;			// event descriptor object

public:
	MTTip(Microtubule*, TrajectoryVector, double, DeterministicQueue*, double);
	~MTTip();

	TipType type();								// returns the tip type (plus or minus)
	double position();							// returns the position of the tip
	double otherPosition();						// returns the position of the other side of the segment
	Segment& segment();							// returns the segment the tip is part of


	void initialize();							// links the tip to the trajectory and determines the first event
	void unlinkFromTrajectory();	//?!?!?! 	// unlinks the tip from the trajectory (called before destruction)
	void switchTrajectory(Trajectory*, Direction, IntersectionItr, bool = true);
												// switches the tip from one trajectory to another


	void locateIntersection();					// locates the next intersection that will be encountered on the trajectory
	void advanceIntersection();					// advances 'nextCollision' by one step
	void determineEvent();						// determines the type of event that takes place at 'nextCollision' and queues it

	void notifyInsert(IntersectionItr&);		// called by trajectory whenever a new intersection is inserted on the trajectory
	void notifyRemove(IntersectionItr&);		// called by trajectory whenever an intersection is removed

private:
	// avoid accidental (expensive) copying of MTTip objects, by declaring private copy constructors without definitions
	MTTip(const MTTip&);
	MTTip& operator=(const MTTip&);

};







class Segment : public DLBaseItem<Segment>
{
public:
	Microtubule* mt;			// pointer to containing microtubule
	Trajectory* const trajectory;	// pointer to containing trajectory
	TrjSegmentTag trajectoryTag;	// Token that links the segment to the trajectory
	double nucleationTime;			// Time of initial nucleation

	Direction dir;					// -+ direction with respect to trajectory orientation

	double start;					// position of minus end
	double end;						// position of plus end
	IntersectionItr startItr;		// iterator pointing to the intersection at or just beyond (for tips) the minus end
	IntersectionItr endItr;			// iterator pointing to the intersection at or just beyond (for tips) the plus end

	Segment(Microtubule*, TrajectoryVector&);
			// constructs a segment as part of a microtubule at a specified vector location
	~Segment();

	double length() {return abs(end-start); }
	bool isLastInMT();
	bool isFirstInMT();
	bool crossesIntersection(IntersectionItr& is) ;
};






class Microtubule : public DLBaseItem<Microtubule>
{
public:
	System* const system;			// pointer to overall system
	double nucleationTime;	// time of initial nucleation event
	int previousUpdateTag;			// time tag at which the positions and lengths were last updated

	MTType type;					// type of microtubule (growing, shrinking)
	MTTip plus;						// the plus tip object, with associated events
	MTTip minus;					// the minus tip object, with associated events
	DLList<Segment>	segments;		// list of segments
	EventDescriptor disappearEvent;	// event descriptor for events that are not directly associated to a single tip (only disappearance)
	
public:
	Microtubule(System*, TrajectoryVector, bool = true);		// nucleates a new microtubule at a specified trajectory position
	~Microtubule();
	bool integrityCheck();						// performs an integrity check on the microtubule

	// status functions
	double length();
	void updateLength(bool forceUpdate = false);			// updates the length to the current system time tag

	// event functions below
	void setDisappearEvent();					// calculates the next 'disappear' event, and queues it
	void handleEvent(const EventDescriptor*);			// gets called to execute an event associated with this MT
	void catastrophe();						// performs a catastrophe
	void rescue();							// performs a rescue
	void harakiri();						// self-destruct event (plus and minus ends meet)
	void wall();							// called when a microtubule hits a region boundary
	void collision();						// executes a plus end collision
	void zipper(Direction dir);					// performs zippering (called from (collision()))
	void crossover();						// called whenever a microtubule crosses an intersection
	void backtrack(MTTip*);						// called when a microtubule retreats across an intersection
	void endOfSegment(MTTip*);					// called when the end of a connected segment is reached (shrinking around a corner/ across a boundary)
	void sever(Segment*,double);						// severs the microtubule at the requested segment and position
	void severAtCross(IntersectionItr is, Segment* cutSeg);		// severs the microtubule at the specified intersection.
	void splitSegmentAtTrajPos(double cutPos, Segment* cutSeg);	// performs the actual splitting of sever... and schedules the next event etc.
	void translatePositionMT2Segment(double& cutPos, Segment*& cutSeg);	 //translates random position at MT to corresponding position at Segment

private:
	// avoid accidental (disastrous) copying of Microtubule objects, by declaring private copy constructors without definitions
	Microtubule(const Microtubule&);
	Microtubule& operator=(const Microtubule&);
};



/******************************* System class *****************************/


class Measurement {
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

	double G_effAdjust_normal;
	double G_effAdjust_special;
};


class ExtensibleHistogram {
	int bins;
	double range;
	int totalCount;
	int* histogram;
public:
	ExtensibleHistogram(int, double = 1.0);
	~ExtensibleHistogram();
	void insert(double);
	void flush(ostream&);
	void initialize(int, double = 1.0);
	bool empty() { return (totalCount == 0); }

private:
	void extendRange();
};

class SingleHistogram {
	ExtensibleHistogram h;
	ofstream file;
	int bins;
	System* system;
public:
	bool valid() { return (bins != 0); }
	void insert(double);
	void save();
	SingleHistogram(System*);
	void setBins(int);
	void setFileName(string);
	void closeFile();
	~SingleHistogram();
};

class MultiAngleHistogram {
	int angles;
	int bins;
	ofstream file;
	System* system;
	ExtensibleHistogram** h;
public:
	bool valid() { return ((angles != 0) && (bins != 0)); }
	void insert(double,double);
	void save();
	MultiAngleHistogram(System*);
	void setBins(int,int);
	void setFileName(string);
	void closeFile();
	~MultiAngleHistogram();

};

class Parameters
/*
 * This object contains all simulation parameters and the functions needed to
 * import/export them from/to the filesystem. Also, parameter verification is
 * performed here.
 * 
 */
{
public:
	System* system;

	// biological parameters
	double vPlus;			// plus end growing velocity 
	double vMin;			// plus end shrinking velocity
	double vTM;				// minus end shrinking velocity
	double kSev;			// severing rate [1/(sec*mu)]
	double kCross;			// severing rate at (occupied) intersections [1/#OccupiedIntersections/sec]
	double kCat;			// spontaneous catastrophe rate [1/sec]
	double kRes;			// spontaneous rescue rate
	double catastropheMultiplier;

	int treadmillingEnabled;
	int severingEnabled;
	int crossSeveringEnabled;
	int crossSeveringTop;
	double crossSeveringStartAngle;

	int forbiddenZones;

	int edgeCatastropheEnabled;
	int edgeCatastropheSmooth;
	double pCatRegularEdge;
	double pCatSpecialEdge;

	// nucleation stuff
	double kNuc;			// nucleation rate [1/(sec*mu^2)]
	NucleationType nucleationType;
	double nucleationAlpha;
	double * nucleationAngles;
  double preSeededSeedDensity;
  double preSeededRate;
  NucleationType preSeededType;
  double preSeededAlpha;

	// interaction stuff

	int zipperingEnabled;				// 1 if zippering is enabled
	int catastrophesEnabled;			// 1 if induced catastrophes are enabled
	int proportionalCatastrophes;		// catastrophe fraction is proportional to collision angle
	double inducedCatastropheFraction;	// catastrophe fraction at perpendicular collision
	double zipFraction;			// fraction of zippering events for th < magicAngle (zipFirst interaction type only)
	double catStartAngle;					// no induced catastrophes below this angle
	double magicAngle;					// switching angle between zippering and crossover
	double c0Value;
	double z0Value;

	InteractionType interactionType;
	BundleType bundleType;
	int discreteAngleNumber;

	// chanLloyd nucleation no longer maintained
	//double chanLloydForwardFraction;		// Fraction nucleated along the existing MT, in the same direction
	//double chanLloydBackwardFraction;		// ... in the reverse direction
	//double chanLloydSidewaysAngle;			// Average angle for sideways nucleations
	//double chanLloydSidewaysWidthToAverageAngle; 	// Half width of the "triangle" (min to avg) or (avg to max). total width is 2x this.
	double nucleationHalfIsotropicDensity;		// Competition between microtubule based and (background) isotropic nucleation (x/(x+HIL))
        int ellipseReducedFreeRate;                  // If true, a fraction of free nucleation events is rejected
        double ellipseReducedFreeRateAcceptFraction;    // This fraction of free nucleation events is accepted. No point in making this >1. 
	double ellipseEpsilon;				// Eccentricity of ellipse
	int ellipseForwardAlongMT;			// Forward/backward nucleations not elliptic, but exactly along MT
	double ellipseLeftFraction;			// Fraction of nucleation events rotated left (by sidewaysAngle degrees)
	double ellipseRightFraction;			// ... right
	double ellipseBackwardFraction;			// ... in the reverse direction (180 degrees)
	double ellipseSidewaysAngle;			// Average angle for sideways nucleations

	
	// geometry stuff
	GeometryType geometry;		// defines the type of geometry (periodic, cylinder, sphere, etc.)
	double geomParam1; 			// 3 parameters describing the size of the geometry (3 ought to be enough)
	double geomParam2;
	double geomParam3;


	// simulation parameters
	unsigned long seed;
	double stopTime;			// end time of the simulation
	double measurementInterval;	// interval between subsequent 'status' events
	double densityLimit;		// microtubule density limit [per micrometer]
	double wallClockLimit;		// run time limit [seconds]
	double memoryLimit;			// memory usage limit [MB]



	// output parameters
	string outputDir;
	int createSubdir;

	double newParameterReadInterval;
	string newParameterFile;

	int movieEnabled;
	double movieFrameInterval;
	int angleHistogramBins;
	int angleHistogramOptical;

	int hiresLengthHistogramBins;
	int loresLengthHistogramBins;
	int loresAngleHistogramBins;
	int loresLifetimeHistogramBins;
	int hiresLifetimeHistogramBins;
	int histogramAverageSamples;


	// theory parameters
	double c0calc;
	double x0calc;
	double z0calc;
	double G_0;
	double adjustedG;
	double lengthFactor;		// l = lengthFactor * L 

public:
	Parameters(System*);			// initialize parameters by reading the file (if not NULL)
	void initialize(const char *);
	bool reinitialize(const char *);
	bool readFromFile(const char*, bool);	// reads parameters from file and checks them
	bool writeToFile();			// writes all parameters to file
	void verifyParameters();	// performes a consistency check on the parameters and adjusts where necessary
	bool calcTheoryParameters();	
	// path outputDir();
};



class System
/* 
 * Master object defining an interacting mt system.
 * 
 */
{
public:
	Geometry* geometry;												// pointer to containing geometry
	DLList<Microtubule,MICROTUBULE_GRANULARITY> growing_mts;		// list of growing MTs
	DLList<Microtubule,MICROTUBULE_GRANULARITY> shrinking_mts;		// list of shrinking MTs
	CompactList<OccupiedIntersection,OCCUPIED_INTERSECTION_GRANULARITY> OccupiedIntersectionList;
	int growingTipsSpecial;
	int growingTipsNormal;

	Parameters p;				// parameter object
	MTRand randomGen;			// random generator

	// global system variables
	double systemTime;			// current system time with respect to offset
	double systemTimeOffset;	// current system time offset
	int currentTimeTag;			// time tag associated with the system time
	double totalLength;			// total length of all microtubules in the system
	int seedsLeft ;				// number of remaining pre-seeded nucleation sites
	int memUsage;				// total memory usage in MB
	bool stopSignal;			// flag that can be set to stop the simulation at the next event
	bool stopping;			// indicates that a stop threshold has been exceeded. No further threshold
								// checks will be performed.

	// simulation statistics
	time_t wallClockStartTime;
	int totalSEventCount;		// cumulative number of stochastic events
	int totalValidDEventCount;	// cumulative number of valid deterministic events
	int totalInvalidDEventCount;// cumulative number of invalid deterministic events
	int totalLengthSeveringCount; // cumulative number of severing events at random positions (sever)
	int totalIntersectionSeveringCount; // cumulative number of severing events at occupied intersections (severAtCross)
	int boundaryCrossingCount;
	int totalZipperCount;
	int totalCrossoverCount;
	int totalInducedCatastropheCount;

	// for memory usage tracking
	int countSegments;
	int countTrajectories;
	int countIntersections;
	int estimateMemoryFootprint();

	// event queues (stochastic and deterministic)
	double nextStochasticEventTime;				// time of next stochastic event
	StochasticEventType nextStochasticEventType;// type of next stochastic event
	DeterministicQueue timeQueue;				// queue object for time-defined events
	DeterministicQueue vPlusQueue;				// queue object for vPlus-defined events 
	double identity(double i) {return i;}		// identity function function
	double vPlusToTime(double);					// distance to time conversion for growing tips
	double timeToVPlus(double);					// time to distance conversion for growing tips
	
	// global event occurances
	double nextStatusEventTime;
	double nextSnapshotEventTime;
	double nextParameterEventTime;

	// manage event descriptor objects
	EventDescriptorIndex EventDescriptorID;							// current tag for event descriptor
	map<EventDescriptorIndex, EventDescriptor*> EventDescriptorMap;	// structure that maps EventDescriptorIDs onto EventDescriptor objects
	EventDescriptorIndex registerEventDescriptor(EventDescriptor*);	// returns a new event descriptor tag
	void unregisterEventDescriptor(EventDescriptorIndex);			// invalidates an event descriptor tag
	EventDescriptor* getEventDescriptor(EventDescriptorIndex);		// returns a pointer to an EventDescriptor through its tag

	// manage event tracking tags
	EventTrackingTag eventID;										// the current event tracking tag
	EventTrackingTag getEventTag() { return eventID++; }			// returns a new event tracking tag

public:
	// initialization/destruction functions
	System(char*);							// initializes the system from a parameter file
	~System();
	bool integrityCheck();					// checks the internal consistency of the variables
	void flushAndReload(bool = true);

	void advanceTime(double nextTime);		// advances the system time to 'nextTime'
	void updateAll(bool forceUpdate = false);						// updates all positions in the system to the system time
	void emergencyBreak();					// writes stored data to file and exits ASAP

	// running and handling events
	void run(double);						// runs the system for a certain amount of simulation time
	void nextEvent(void);					// gets the first event from the queue and handles it
	void handleGlobalEvent(DeterministicEvent&);
	void determineStochasticEvent();		// function that determines and queues the next stochastic event
	void handleNucleationEvent(bool=false);
	void handleSeveringEvent();
	void handleSeveringAtCrossEvent();
	void handleRescueEvent();
	void handleCatastropheEvent();
	void randomPositionOnMicrotubule(double& cutLength, Segment*& cutSeg);

	CollisionType collisionResult(double, int, int);	// returns the result of a collision at a specified angle, and bundle occupancies
	void collisionProbabilities(double, double&, double&);
		// returns the probabilities for induced catastrophes and zippering, as a function of the angle
	// binomial table. Necessary for bundle collision type: multicollisions
	double binomialTable[MAXBINOM][MAXBINOM];
	void makeBinomialTable(void);
	double multiPcross(int Npar, int Ncoll, double xSingle, double zSingle);
	double multiPzip(int Npar, int Ncoll, double xSingle, double zSingle);




	ofstream parameterFile;

	ofstream movieFile;
	ofstream measurementFile;
	ofstream angleHistogramFile;
	ofstream angleHistogramOpticalFile;
	ofstream angleLengthFile;
	ofstream angleNumberFile;

	deque<Measurement> measurementHistory;
	deque<double *> angleHistory;
	deque<double *> angleHistoryOptical;
	deque<double *> angleLengthHistory;
	deque<int *> angleNumberHistory;

	SingleHistogram mtLengthHistogram;
	SingleHistogram mtLifetimeHistogram;
	MultiAngleHistogram segAngleLengthHistogram;
	MultiAngleHistogram segAngleLifetimeHistogram;
	int histogramAverageCount;

	bool dataSaved;
	void initializeOutput();
	void performMeasurement();
	void writeMeasurementsToFile(int = 0);
	void closeFiles();

	void removeOccupiedIntersection(Intersection& is);	// removes pointers from is and its mirror
	void addOccupiedIntersection(IntersectionItr is);   	// creates pointers for is and its mirror 
								// OccupiedIntersection contains only one pointer: to the one "on top" (to be cut).

private:
	// avoid accidental (expensive) copying of Trajectory objects, by declaring private copy constructors without definitions
	System(const System&);
	System& operator=(const System&);
};


ostream& operator<<(ostream&, const Measurement);
void writeMeasurementDescriptors(ostream&);

/********************* Geometry specializations *******************************/

class Periodic : public Geometry
/*
 * Periodic boundary conditions on a rectangular surface
 * 
 */
{
	const Coord2D size;				// x and y dimensions
public:
	Periodic(Coord2D, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[], bool=false);
	void outputSnapshot(ostream&);
	
	
	// this function is specific to the periodic geometry
	void calculateDiscreteAngleLengths(double[], int[]);
};

class Grid : public Geometry
/*
 * Periodic boundary conditions on a rectangular surface
 * 
 */
{
	const Coord2D size;				// x and y dimensions
	const int xNumber;
	int yNumber;
public:
	Grid(Coord2D, int, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	
	// this function is specific to the periodic geometry
	void calculateDiscreteAngleLengths(double[], int[]);
};


class Box : public Geometry
/*
 * Rectangular 3D box geometry
 */
{
	const double x;
	const double y;
	const double z;
public:
	Box(double, double, double, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
};

class Cylinder : public Geometry
/*
 * 3D cylinder geometry
 */
{
	const double radius;
	const double length;
public:
	Cylinder(double, double, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
};

class GridCylinder : public Geometry
/*
 * Cylinder with grid on the main surface
 * 
 */
{
	const double radius;
	const double length;
	const int number;
	int gridRNumber;
	double gridRSize;
	double gridLenSize;
public:
	GridCylinder(double, double, int, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	
};

class Pancake : public Geometry
/*
 * 3D cylinder geometry
 */
{
	const double radius;
public:
	Pancake(double, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
};

class Wormhole : public Geometry
/* 
 * Wormhole geometry. This is an adaptation of the rectangle with periodic boundary conditions 
 * that approximates one aspect of the spherocylinder: the end caps. Microtubules that exit the 
 * rectangle on the left or right (x=0 or x=max) come back at the same side, under a flipped
 * angle and translated by half a box size.
 * 
 * 
 */
{
	const Coord2D size;
public:
	Wormhole(Coord2D, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);

};

/********************************* Region specializations ***************************/


class Cartesian : public Region
{
public:
	Cartesian(RegionType, Geometry*, int, double, RegionAccountingType, RegionParameterType);

	// virtual function implementations	
	virtual SurfaceVector randomSurfaceVector() = 0;				// produces a random surface vector inside the region
	virtual void getTrajectoryCoordinates(SurfaceVector&, double&, TrajectoryVector&) = 0;
														// Used for nucleation:
														// adjusts the surface vector to become the base for a trajectoryvector, and 
														// returns the total length and trajectory vector (minus trajectory label).
	void makeIntersectionList(Trajectory*);				// calculates the intersections of the given trajectory with all others in the region
	void translateVector(SurfaceVector&, const double);	// translates the surface vector by a given distance
	double intersectionAngle(Trajectory*, Trajectory*);	// determines the intersection angle between two trajectories
	void getOrderParameters(OrderParameters&);
	void getOrderParametersRawFlat(OrderParametersRaw&, double[3][2]);
	void getOrderParametersRawCylinder(OrderParametersRaw&, double, double, double);
	double calculateHistogram(double[], bool optical = false, bool partOfSequence = false); 	// return value: localLength. 
	void outputSnapshotOffset(ostream&, double, double);
	void outputSnapshot(ostream&);
};

class GRectangle : public Cartesian
/*
 * Simple rectangle.
 * 
 */
{
	const Coord2D size;				// x and y sizes
public:
	GRectangle(Coord2D, Geometry*, int, RegionAccountingType = r_accounting_normal, RegionParameterType  = r_param_normal);

	// virtual function implementations	
	SurfaceVector randomSurfaceVector();				// produces a random surface vector inside the region
	void getTrajectoryCoordinates(SurfaceVector&, double&, TrajectoryVector&);
														// Used for nucleation:
														// adjusts the surface vector to become the base for a trajectoryvector, and 
														// returns the total length and trajectory vector (minus trajectory label).

	// Special function for the rectangle object - not required for other regions
	RectangleSide locateSide(double&, double&);			// given an x and y coordinate, returns the closest side of the rectangle
};



class Disc : public Cartesian
/*
 * Simple disc.
 * 
 */
{
	const double radius;				// x and y sizes
public:
	Disc(double, Geometry*, int, RegionAccountingType = r_accounting_normal, RegionParameterType  = r_param_normal);

	// virtual function implementations	
	SurfaceVector randomSurfaceVector();				// produces a random surface vector inside the region
	void getTrajectoryCoordinates(SurfaceVector&, double&, TrajectoryVector&);
														// Used for nucleation:
														// adjusts the surface vector to become the base for a trajectoryvector, and 
														// returns the total length and trajectory vector (minus trajectory label).
};



#ifndef NO_INLINE
#include "inline.cpp"
#endif


#endif /*CORTICALSIM_H_*/
