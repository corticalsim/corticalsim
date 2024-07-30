#ifndef CORTICALSIM_H_
#define CORTICALSIM_H_

/* corticalSim - Cortical microtubule simulator
 * Simon Tindemans and Eva Deinum, FOM Institute AMOLF, 2008-2014
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
 * * If the tubulin pool size is decreased on the fly, care should be taken that the actual density doesn't exceed the
 *   pool size. If this occurs, the program behaviour is undefined.
 * * In microtubule.cpp and system.cpp, compiler warnings are (should be) issued for using the 'this' pointer within
 *   the initializer list, because the 'this' pointer cannot be used until the initialization has completed. The 
 *   current use is ok, because the pointer is only used to store the address of the parent object. NOTE: when making 
 *   changes to the code, take care NOT to use this pointer within the constructors of other objects!

 TO DO:
 - remove question marks and exclamations in intersection ordering code (made by eva?)
 */


#define PROGRAM_VERSION "2.0 beta"

//#define VAR_CAT // if defined, it is possible to have variable r_cat by sine function
//#define BAND_CAT // if defined, it is possible to have variable r_cat by distinct bands   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

#define COMPILER_RECURSIVE_SUPPORT    // if defined, defines 'class Intersection' using an incomplete type.
                                        // This leads to a large reduction in redundant memory allocations for
                                        // compilers that support it (gcc, VC++, but *not* CLANG/Xcode)
                                        // DBG_ASSERT code fails without this support.
//#define NO_INLINE				// if defined, avoids inlining the functions in inline.cpp. This 
								// carries a performance penalty, but may aid debugging
//#define NO_INTEGRITY_CHECK	// if defined, disables the system integrity checks performed
								// at every measurement. This may speed up execution if you perform 
								// measurements at a high rate, but makes the code less aware of
								// potential bugs. Not recommended.
//#define SLEEZY // if defined, integrity check may be carried out, but code continues running upon failure. Not recommended.

//#define MTBASED_NUCLEATION_PROBABILITY               //if defined, calculates the probability of MT based nucleation as a function of the local density of tubulin

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
//#define DBG_SPATIAL
//#define DBG_VELOCITY
//#define DBG_ASTER
	
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
#include <boost/container/map.hpp>
#include <boost/variant.hpp>
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
const int METATRAJECTORY_GRANULARITY = 512; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// TODO CLARIFY THIS!!!!!!!
const int INTERSECTION_GRANULARITY = 4096;
const int EVENT_GRANULARITY = 1024;
const int OCCUPIED_INTERSECTION_GRANULARITY = 256;

const int MAX_HISTORY_SIZE = 256; // maximum size of measurement cache (increase to minimize disk access)
const int MIN_HISTORY_SIZE = 64;  // minimum size of measurement cache (minimum may be needed for trending) 

// simulation parameters
const int POSITION_CACHE_SIZE = 32768;
const int CLOCK_POLLING_INTERVAL = 10000;
const int MEMORY_POLLING_INTERVAL = 200;
const int QUEUE_FLUSH_INTERVAL = 1000000;//1000000;			//note that this is related to the total number of events...

const int NUCLEATION_DISCRETIZATION_STEPS = 512;



typedef enum {g_periodic, g_grid, g_wormhole, g_box, g_cylinder, g_gridcylinder, g_spherocylinder, g_sphere, g_pancake, g_COUNT_LAST} GeometryType;
extern string GeometryTypeText[]; // = {"periodic", "grid", "wormhole", "box", "cylinder", "gridcylinder", "spherocylinder", "sphere", "pancake"};

typedef enum {r_rectangle, r_disc, r_dome} RegionType;
extern string RegionTypeText[]; // = {"rectangle", "disc", "dome"};

typedef enum {o_x, o_y, o_z, o_COUNT_LAST} OrientationType;
extern string OrientationTypeText[]; // = {"x", "y", "z"};

typedef enum {sh_none, sh_x, sh_y,  sh_xy, sh_COUNT_LAST} SpatialHistogramType;
extern string SpatialHistogramTypeText[]; // = {"none", "x", "y", "xy" };


typedef enum {forward=1, backward=-1} Direction;
typedef enum {s_top, s_left, s_bottom, s_right} RectangleSide;

typedef enum {ev_none, ev_wall=1, ev_collision, ev_backtrack, ev_end_of_segment, ev_disappear, ev_deflection} DeterministicEventType;
typedef enum {status=1, snapshot, stop, parameter_change} GlobalEventType;
typedef enum {catastrophe=1, rescue, katanin, severingAtCross, nucleation, preSeededNucleation,extraRescue} StochasticEventType;
typedef enum {ct_zipper=1, ct_crossover, ct_inducedCatastrophe} CollisionType;

//typedef enum {nuc_isotropic, nuc_aster, nuc_biased, nuc_discreteAngles, nuc_chanLloyd, nuc_chanLloydRandomPosition, nuc_chanLloydIsotropic, nuc_ellipse,  nuc_COUNT_LAST} NucleationType;	
//extern string NucleationTypeText[]; // = {"isotropic", "biased", "discreteAngles", "chanLloyd", "chanLloydRandomPosition", "chanLloydIsotropic", "ellipse" };
typedef enum {nuc_isotropic, nuc_aster, nuc_biased, nuc_discreteAngles, nuc_ellipse, nuc_COUNT_LAST} NucleationType;	
typedef enum {nbias_cross, nbias_ellipsePolar, nbias_ellipseApolar, nbias_COUNT_LAST} NucleationBiasType;	
extern string NucleationTypeText[]; // = {"isotropic", "biased", "discreteAngles", "ellipse" };
extern string NucleationBiasTypeText[]; // = {"cross", "ellipsePolar", "ellipseApolar" };

typedef enum {int_zipFirst,  int_catFirst, int_minimalFourier, int_COUNT_LAST} InteractionType;
extern string InteractionTypeText[];// = {"zipFirst", "catFirst", "minimalFourier" };

typedef enum {bdl_simple,  bdl_sticky, bdl_noZip, bdl_multiCollision, bdl_Ncollision, bdl_COUNT_LAST} BundleType;
extern string BundleTypeText[];// = {"simple", "sticky", "noZip", "multiCollision", "Ncollision" };

typedef enum {t_minus, t_plus} TipType;
typedef enum {mt_growing, mt_shrinking} MTType;

typedef enum {graphics_mt, graphics_trajectory, graphics_plus, graphics_minus} GraphicsType;

typedef enum { r_accounting_normal, r_accounting_dontcount } RegionAccountingType;
typedef enum { r_param_normal, r_param_modified } RegionParameterType;

typedef enum { p_bar, p_cross, p_band, p_singleBand, p_COUNT_LAST } PatternType;
extern string PatternTypeText[];

typedef enum { r_cos, r_COUNT_LAST } extraRescueType;

extern string extraRescueTypeText[];

class System;
class Microtubule;
class Segment;
class MTTip;
class DeterministicEvent;
class Parameters;

class Geometry;
class Region;
class Trajectory;
class MetaTrajectory;
class Intersection;
class OccupiedIntersection;
class OneDSpatialMeasurement;


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



struct Aster 
{       /* Very similar concept to SurfaceVector, except it has nor distinct angles.
	 * 
	 * Note that the names x, y and angle are chosen to reflect their meaning for 
	 * for 2D cartesian surfaces, but their meaning may differ depending on the region type
	 */
  	int nor;
  	double x;
  	double y;
  	vector<double> angles;
  	Region* region;
};


struct NucSpotCandidate
{       /* 
	 * Defines a possible nucleation spot at a certain SurfaceVector
	 */
	Trajectory* motherTrajectory;
	double posOnMother;
	int numOfSegmentsInBundle;
	SurfaceVector position;
	double distanceFromAster;
};


struct localDensity
{	/*
	 * Defines the local density, nearby an insertion point
	*/
	vector<Region*> nearbyRegions;                                  // list of regions close to the insertion point to calculate the local density /////////////////////////////////////////////////////////////
	//double nearbyTotalLength;
	double nearbyDensity;
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
  double localL;
  double localLOpt;

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




/* NOTE:
 The definition of class Intersection contains an IntersectionItr, which is an iterator to a multimap that contains Intersection itself. According to the C++ standard, the outcome is not defined and some compilers (notably Xcode) reject it. This is remedied using the boost::recursive_wrapper, but it results in large numbers of memory allocations. If the compiler supports multimap iterator definitions using incomplete types, defining COMPILER_RECURSIVE_SUPPORT avoids this overhead.
 */
#ifdef COMPILER_RECURSIVE_SUPPORT
typedef multimap<double, Intersection, std::less<double>, boost::fast_pool_allocator<std::pair<double, Intersection> > > IntersectionMap;
#define ISREF second
#else
typedef multimap<double, boost::recursive_wrapper<Intersection>, std::less<double>, boost::fast_pool_allocator<std::pair<double, boost::recursive_wrapper<Intersection> > > > IntersectionMap;
#define ISREF second.get()
#endif
typedef IntersectionMap::iterator IntersectionItr;


class Intersection
/*
 * The intersection object describes the intersection of one trajectory with another, from one side.
 * Since it is created twice for every intersection, it is by far the most memory-intensive object
 * and should be kept as small as possible.
 */
{
public:
    Intersection() : occupancy(0), otherTrajectory(nullptr), occupiedListPtr(nullptr) {};
    
    int occupancy;					// number of MTs that are crossing the intersection along this trajectory
	Trajectory* otherTrajectory;	// pointer to the other trajectory
	IntersectionItr mirror;			// pointer to the mirroring Intersection structure of the other trajectory
	// TODO: maybe add corner count
	OccupiedIntersection* occupiedListPtr;
};


class OccupiedIntersection : public CompactListItem<OccupiedIntersection>
{
public:
    OccupiedIntersection(IntersectionItr is) : intersectionToCut(is) {};
    IntersectionItr intersectionToCut;                                    
    //IntersectionItr metaIntersection;		                                  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};

class Trajectory : public DLBaseItem<Trajectory>
/*
 * The trajectory is the basic geometrical object. Tips and segments associate with a trajectory
 * and move/lie alongside it. Trajectory intersections determine the collision points.
 * 
 */
{
	friend class Region;
	friend class MetaTrajectory;

public:
	const SurfaceVector		base;			// coordinates and angle of the base of the trajectory
	const double 			length;			// total length of the trajectory, from one boundary to another
	TrajectoryVector		prevTr;			// continuation of the Trajectory at the base
	TrajectoryVector		nextTr;			// continuation of the Trajectory at the far end

	double					prevTrCosAngle;	// Cosine of the 3D angle with the previous trajectory (for edge catastrophes)
	double					nextTrCosAngle; // Cosine of the 3D angle with the next trajectory (for edge catastrophes)

    IntersectionMap intersections;
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
	//void coveredLengthSplit(double * length, double* Olength, double* cuts, int n, int); // combination of coveredLength ("Olength") and segmentLength ("length"), divided in n pieces by cuts
	void coveredLengthSplit(double * length, double* Olength, double* count,double* cuts, int n, int); // combination of coveredLength ("Olength") and segmentLength ("length"), divided in n pieces by cuts // count doesn't work yet.
  int numberOfSegmentsPassingPoint(double position); // find number of segments in trajectory that pass a certain point, excluding segments that precisely reach it.
  int numberOfOtherSegmentsPassingIntersection(IntersectionItr& refis, Segment* refseg); // Find number of segments passing an intersection excluding the reference segment
  int numberOfSegmentsDeflectingAtIntersection(IntersectionItr& is, Segment* refseg); // get all segments that switch trajectory at intersection

private:
	// avoid accidental (expensive) copying of Trajectory objects, by declaring private copy constructors without definitions
	Trajectory(const Trajectory&);
	Trajectory& operator=(const Trajectory&);
};



class MetaTrajectory : public DLBaseItem<MetaTrajectory>
{
	//friend class Region;
	//friend class Geometry;

public:
	System* const system;                           //pointer to the overall system: I need it cause the meta trajectory needs to know everything that is going on
	const double maxDistance;                          //max distance before assuming that the closest MT is too far for MT-based nucleation TODO: add it to the parameter list 
	const SurfaceVector startOfMetaTrajectory;           //location of the centre of the aster and angle
	bool doesItIntersectOccupied;
	Direction firstDirection;
	#ifdef MTBASED_NUCLEATION_PROBABILITY
	vector<Region*> crossedRegions; 
	#endif
	
	bool getFindMT() { return doesItIntersectOccupied; }                    //returns whether meta trajectory meets a MT lattice  
	vector<Trajectory*> trajectoriesOfRay;         //pointers to all trajectories that form the meta trajectory
	
	NucSpotCandidate nsc;
	double metaTrajectoryLength;                   //total length of the meta trajectory
		
	MetaTrajectory(System* const, double /*maxDist*/, SurfaceVector);              //any SurfaceVector, but it should take the SurfaceVector from the Aster
	~MetaTrajectory();
	
	void addFirstTrajectory();                           //add the first trajectory to the vector
	void removeTrajectories();                                //remove all trajectories of the Meta Trajectory
	
};





class Region
/* 
 * Virtual base class for a geometry patch with a 2D coordinate system. Within a region,
 * trajectories may only intersect once.
 * 
 * Region functions are only called by geometry or (meta)trajectory member functions.
 */
{
	friend class Geometry;
	friend class Trajectory;
	friend class MetaTrajectory;

public:
	Geometry* const 	geometry;			// pointer to containing geometry
	const int geometryRegionIndex;			// region index corresponding to this region within the geometry.
	const RegionType	type;				// type of region (rectangle, disc, hemisphere)
	const double 		area;				// area of the region
	const RegionAccountingType accountingType;	// whether segments in this region are part of the 'main' region or not
	const RegionParameterType parameterType;

	DLList<Trajectory, TRAJECTORY_GRANULARITY> trajectories;					// list of trajectories on the region
	

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
	virtual Aster createAster(int, SurfaceVector) =0;                       // generates random angles for the aster

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
	friend class MetaTrajectory;

public:	
	const GeometryType 	type;					// the geometry type (periodic, sphere, etc.)
	const double 		area;					// total area of all surfaces in the geometry
	System* const		system;					// pointer to the containing system
//protected:
	vector<Region*> regions;					// list of the Region-derived objects that make up the geometry
	
	DLList<MetaTrajectory, METATRAJECTORY_GRANULARITY> metaTrajectories;

public:
	Geometry(GeometryType g, System* s, double a) : type(g), area(a), system(s) {};
	// virtual destructor implies virtual class
	virtual ~Geometry(){ return; };
	bool integrityCheck();						// performs an integrity check of the trajectory network
	double opticalLength();
	int trajectoryCount();
	virtual void getOrderParameters(OrderParameters&) = 0;
  virtual void getLocalOrderParameters(vector<OrderParameters>&) = 0;
	virtual double calculateHistogram(double[],bool=false) = 0;
	virtual void outputSnapshot(ostream&) = 0;
	virtual void getOneDmeasurement(OneDSpatialMeasurement&,double) = 0;
  virtual void fixRegion(SurfaceVector&) = 0; // Used for nucleation if x,y position is shifted independent from trajectory
                            // Checks if x,y falls outside region and adjusts region and x,y accordingly if needed.
                            // Implemented for grid and gridcylinder only. 
  virtual void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double) = 0; // Implemented for grid and gridcylinder only. 


public:
	Aster createAster(int, SurfaceVector);
	void removeMetaTrajectory(MetaTrajectory*);
	SurfaceVector randomSurfaceVector();		// returns a random surface vector
	TrajectoryVector createTrajectory(const SurfaceVector&, bool scanExisting = false);
												// creates a trajectory at a given surface vector and returns a trajectory vector
	TrajectoryVector createAndLinkTrajectory(const SurfaceVector&, Trajectory*, Direction, double);
	

  // Turn xpos on gridcell into xpos on domain
  virtual double xPosGridToDomain(double xPos, int ridx) {return xPos;}
  // Turn ypos on gridcell into xpos on domain
  virtual double yPosGridToDomain(double yPos, int ridx) {return yPos;}

  // Get total area used in spatial histograms (no caps)
  virtual double getAreaForSpatialHist(double geomPar1, double geomPar2) {return geomPar1*geomPar2;}

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
	friend class MetaTrajectory;
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
	bool crossesPoint(double position) ;
	bool deflectsAtIntersection(IntersectionItr& is, Direction refdir) ;
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
  void handleDeflection(); // perform deflection

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
  vector<OrderParameters> local_order;

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
  int nucleationCount;
  int deflectionCount;

	double G_effAdjust_normal;
	double G_effAdjust_special;
	double G_effMeasured;

	double G_effAdjust_min;
	double G_effAdjust_max;
	double G_effAdjust_S2angle;

#ifdef BAND_CAT
  double G_effAdjust_band;
  double G_effAdjust_gap;
#endif

  double nOccupiedNCs;

};

class OneDSpatialMeasurement {
private:
	ofstream file;
	double binSurface;
	double writeRaw;
public:
	int bins;
	double * cuts;
	OrientationType ori ;
	bool includeCaps;
	OrderParametersRaw * opHist;
	double * nHist;
	OneDSpatialMeasurement (int,bool,OrientationType,double,double,bool);
	~OneDSpatialMeasurement();
	void process();
	void writeToFile(double);
	void writeHeader();
	void writeLine(double,int);
	void setFileName(string);
	void closeFile();
	void reset();
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
	double vPlus;			// plus end growing velocity (for infinite tubulin pool)
	double vMin;			// plus end shrinking velocity
	double vTM;				// minus end shrinking velocity
	double kSev;			// severing rate [1/(sec*mu)]
	double kCross;			// severing rate at (occupied) intersections [1/#OccupiedIntersections/sec]
#ifdef VAR_CAT
	OrientationType kCatOrient;			// if spontaneous catastrophe rate varies:  along x / y / z axis (z not implemented)
	double kCatBands;			// if spontaneous catastrophe rate varies along y-axis: # periods of sin function
	double kCatAlpha;			// if spontaneous catastrophe rate varies along y-axis: amplitude
	double invKCatPeriod;			// if spontaneous catastrophe rate varies along y-axis: admin parameter
  double catMin;
  double catMax;
  double Gmin;
  double Gmax;
#endif
double Gmax;                 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef BAND_CAT
	OrientationType kCatOrient;			// if spontaneous catastrophe rate varies:  along x / y / z axis (z not implemented)
  double catMin;
  double catMax;
  double Gmin;
  double Gmax;
  double kCat[2]; // for zone based bands; ultimately with "all" parameters zone dependent: [0] = band ; [1] = gap
  double bandGapWidth[2];  // for zone based bands; ultimately with "all" parameters zone dependent: [0] = band ; [1] = gap
  int nSpirals; // if spirals along cylinder (orientation = x), then n continuous spirals (0 -> rings)
  double spiralPitch; // internal parameter based on nSpirals, orientation and geometry dimensions
  bool calcNucBiasAngle; // if true: calculate nucleationBiasAngle such that it follows the spiral path. 
  double wrapLength; // internal constant for spirals: length of domain to wrap (or 1/meaningless)
  double projectedPeriod; // projection along x/y/z axis of bandGapWidth[0] + [1] 
  double projectedBand; // projection along x/y/z axis of bandGapWidth[0] 
  bool redistributeNucleationOverBands;
  double redistributeShift;
  double redistributeShiftMax;
  int redistributeShiftNumber;
  PatternType catPattern;
  bool reducedGapNucleation; // switch: reduce nucleation if chosen inside gap
  double gapNucleationAcceptFraction; // if reducedGapNucleation, kNuc_gap = kNuc[applicable] * gapNucleationAcceptFraction
#else
	double kCat;			// spontaneous catastrophe rate [1/sec]
#endif
	double kRes;			// spontaneous rescue rate
	double kResExtraMax;			// extra spontaneous rescue rate
	double kResExtraAngle;			// angle of highest rescue rate
	extraRescueType extraRescueFunction;    // function used for angle dependent rescues. Initially implementing only "cos"
	double poolDensity;			// size of tubulin pool [mu]
	double catastropheMultiplier;

	int treadmillingEnabled;
	int severingEnabled;
	int crossSeveringEnabled;
	int crossSeveringTop;
  double crossSeveringTopFraction;
	double crossSeveringStartAngle;
	int restrictedPool;

	int forbiddenZones;

	int edgeCatastropheEnabled;
	int edgeCatastropheSmooth;
	double pCatRegularEdge;
	double pCatSpecialEdge;

	// nucleation stuff
	double kNuc;			// nucleation rate [1/(sec*mu^2)]
	NucleationType nucleationType;
	NucleationBiasType nucleationBiasType;
	double nucleationAlpha;
	double nucleationBiasAngle;
	double * nucleationAngles;
  double preSeededSeedDensity;
  double preSeededRate;
  NucleationType preSeededType;
  double preSeededAlpha;
  bool useRegionalNucleationSaturation;
  double f_unbound;
  bool useDoubleNucleationSaturation;
  double rn_targetratio;
  double occupancytimeNC;
  int nNCmax;
  double rn_base0;
	double gridAspectratio;
  bool useAbsoluteBidirectionalNucBias;
  double absoluteNucBiasAngle;
  double nucBiasVariance;
  
  // nucleation with aster
  int numberOfRays;
  double maxMetaTrajectoryLength;
  double diffusionCoefficient;
  double rvRejectionUnbound;
  double rvRejectionBound;
  double unboundNucleationRate;

  // deflection stuff
  bool useMTdeflection;
  double deflectionStepsize;
  double deflectionMaxAngle;
  double deflectionMinAngle;
  bool deflectionBundleCompensation;
  bool deflectionBundleProtection;
  bool useBundleTracking;
  double bundleTrackMaxAngle;
  double biasDeflectionRight;


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
	double ellipseEpsilonAlongMT;				// Eccentricity of ellipse for (anti)parallel nucleations
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

	SpatialHistogramType spatialHistogramType;
	int spatialHistogramBinsX;
	int spatialHistogramBinsY;
	bool spatialHistogramCountCaps;
	bool spatialHistogramWriteRaw;
	bool regionalOutputQuantities;

  string outputNucPos;


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
  double calcG(double,double); // ugly: uses default parameters or alternative value of kCat / kRes. UPDATE WHEN NEEDED!
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
  int totalNucleationCount;
  int totalDeflectionCount;
  	int totalUnboundNucleationCount;
  	int totalMTbasedNucleationCount;

	// for memory usage tracking
	int countSegments;
	int countTrajectories;
	int countIntersections;
	int estimateMemoryFootprint();

	// event queues (stochastic and deterministic)
	double nextStochasticEventTime;				// time of next stochastic event
	StochasticEventType nextStochasticEventType;// type of next stochastic event
	DeterministicQueue timeQueue;				// queue object for time-defined events
	DeterministicQueue vPlusQueue;				// queue object for vPlus-defined events (regulated by tubulin pool size)
	double identity(double i) {return i;}		// identity function function
	double vPlusToTime(double);					// distance to time conversion for growing tips
	double timeToVPlus(double);					// time to distance conversion for growing tips
	
	// global event occurances
	double nextStatusEventTime;
	double nextSnapshotEventTime;
	double nextParameterEventTime;

  // occupied nucleation time queue
  deque<double> occupiedNCs;

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
  bool handleReducedGapNucleation(SurfaceVector& sv, bool useIsotropic, Segment* nucSeg, double* pospointer);
  void handleBoundNucleation(SurfaceVector& sv, TrajectoryVector& tv, Segment* nucSeg, double& pos);
  double handleBiasedNucleationAngle(SurfaceVector& sv, double overrideAngle, bool preSeeded);
	void handleNucleationEvent(bool=false);
	void handleSeveringEvent();
	void handleSeveringAtCrossEvent();
	void handleRescueEvent(bool=false);
	void handleCatastropheEvent();
	void randomPositionOnMicrotubule(double& cutLength, Segment*& cutSeg);
	void randomSegmentAtMetaIntersection(double posOnTrajectory, Trajectory* tr, int occupancy, Segment*& seg);
  bool handleRegionalNucleationSaturation(int ridx, double& randomPos, Segment*& randomSeg);
  void randomPositionOnMicrotubuleInRegion(int ridx, double cutLength, double& randomPos, Segment*& randomSeg);
  void selectRandomRegionProportionalBy(string quantitytype, double quantitymax, int& ridx, double& cutLength);

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
	#ifdef MTBASED_NUCLEATION_PROBABILITY
	ofstream MTbasedProbDensityFile;
	ofstream unboundProbDensityFile;
	#endif

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
	
	OneDSpatialMeasurement * histX;
	OneDSpatialMeasurement * histY;

  vector<double> nucleationXpositions;
  vector<double> nucleationYpositions;

	bool dataSaved;
	void initializeOutput();
	void performMeasurement();
	void writeMeasurementsToFile(int = 0);
	void closeFiles();

	void removeOccupiedIntersection(Intersection& is);	// removes pointers from is and its mirror
	void addOccupiedIntersection(IntersectionItr is);   	// creates pointers for is and its mirror 
								// OccupiedIntersection contains only one pointer: to the one "on top" (to be cut).

  void initDoubleSatPars(); // initialize extra parameters for double saturating nucleation
  
  	#ifdef MTBASED_NUCLEATION_PROBABILITY
    	vector<double> regionalDensityVector;
    	vector<double> nearbyDensityVector;
    	vector<double> globalDensityVector;
    	vector<double> regionalDensityVectorUnbound;
    	vector<double> nearbyDensityVectorUnbound;
    	vector<double> globalDensityVectorUnbound;
    	#endif

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
public:
	const Coord2D size;				// x and y dimensions
	Periodic(Coord2D, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	void getLocalOrderParameters(vector<OrderParameters>&) {return;};
	double calculateHistogram(double[], bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);
	
	
	// this function is specific to the periodic geometry
	void calculateDiscreteAngleLengths(double[], int[]);
};

class Grid : public Geometry
/*
 * Periodic boundary conditions on a rectangular surface
 * 
 */
{
	friend class MetaTrajectory;
	
public:
	const Coord2D size;				// x and y dimensions
	const int xNumber;
	int yNumber;
	Grid(Coord2D, int, double, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	void getLocalOrderParameters(vector<OrderParameters>&) {return;};
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);
  // Turn xpos on gridcell into xpos on domain
  virtual double xPosGridToDomain(double xPos, int ridx) {
    return xPos + size.x * (-0.5 + 1./xNumber * (ridx%xNumber + 0.5));  
  }
  // Turn ypos on gridcell into xpos on domain
  virtual double yPosGridToDomain(double yPos, int ridx) {
    //return yPos + size.y * (0.5 - 1./yNumber * (ridx/xNumber));
    return yPos + size.y * (0.5 - 1./yNumber * (ridx/xNumber + 0.5));
  }
	
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
	void getLocalOrderParameters(vector<OrderParameters>&) {return;};
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);
};

class Cylinder : public Geometry
/*
 * 3D cylinder geometry
 */
{
public:
	const double radius;
	const double length;
	Cylinder(double, double, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
	void getLocalOrderParameters(vector<OrderParameters>&) {return;};
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);

  // Get total area used in spatial histograms (no caps)
  virtual double getAreaForSpatialHist(double geomPar1, double geomPar2) {return geomPar1*2*geomPar2*PI;}
};

class GridCylinder : public Geometry
/*
 * Cylinder with grid on the main surface
 * 
 */
{
public:
	const double radius;
	const double length;
	double gridLenSize;
	const int number;
	int gridRNumber;
	double gridRSize;
	GridCylinder(double, double, int, double, System*);
	TrajectoryVector extendTrajectory(Trajectory*, Direction);
	void getOrderParameters(OrderParameters&);
  void getLocalOrderParameters(vector<OrderParameters>&);
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);
  // Turn xpos on gridcell into xpos on domain
  virtual double xPosGridToDomain(double xPos, int ridx) {
    return xPos + ((ridx-2)%number-(number-1)/2.)*gridLenSize; 
  }
  // Turn ypos on gridcell into xpos on domain
  virtual double yPosGridToDomain(double yPos, int ridx) {
    return yPos - ((ridx-2)/number-(gridRNumber-1)/2.)*gridRSize;
  }

  // Get total area used in spatial histograms (no caps)
  virtual double getAreaForSpatialHist(double geomPar1, double geomPar2) {return geomPar1*2*geomPar2*PI;}
	
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
	void getLocalOrderParameters(vector<OrderParameters>&) {return;};
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);
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
	void getLocalOrderParameters(vector<OrderParameters>&) {return;};
	double calculateHistogram(double[],bool=false);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
  void fixRegion(SurfaceVector&);
  void shiftSurfaceVectorNoWrap(SurfaceVector&,double,double,double,double);

};

/********************************* Region specializations ***************************/


class Cartesian : public Region
{
public:
	Cartesian(RegionType, Geometry*, int, double, RegionAccountingType, RegionParameterType);

	// virtual function implementations	
	virtual SurfaceVector randomSurfaceVector() = 0;				// produces a random surface vector inside the region
	virtual void getTrajectoryCoordinates(SurfaceVector&, double&, TrajectoryVector&) = 0;
	virtual Aster createAster(int, SurfaceVector) = 0;
														// Used for nucleation:
														// adjusts the surface vector to become the base for a trajectoryvector, and 
														// returns the total length and trajectory vector (minus trajectory label).
	void makeIntersectionList(Trajectory*);				// calculates the intersections of the given trajectory with all others in the region
	void translateVector(SurfaceVector&, const double);	// translates the surface vector by a given distance
	double intersectionAngle(Trajectory*, Trajectory*);	// determines the intersection angle between two trajectories
	void getOrderParameters(OrderParameters&);
	void getOrderParametersRawFlat(OrderParametersRaw&, double[3][2]);
	void getOrderParametersRawFlatBinned(OneDSpatialMeasurement&, double[3][2],  double* , int,int,int,OrientationType,double*);
	void transformCuts(double *cutsIn, int  binsIn, double * cutsOut, const SurfaceVector & sVec, OrientationType mOri );
	void getOrderParametersRawCylinder(OrderParametersRaw&, double, double, double);
	double calculateHistogram(double[], bool optical = false, bool partOfSequence = false); 	// return value: localLength. 
	void outputSnapshotOffset(ostream&, double, double);
	void outputSnapshot(ostream&);
	void getOneDmeasurement(OneDSpatialMeasurement&, double);
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
	Aster createAster(int, SurfaceVector);
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
	Aster createAster(int, SurfaceVector);
	void getTrajectoryCoordinates(SurfaceVector&, double&, TrajectoryVector&);
														// Used for nucleation:
														// adjusts the surface vector to become the base for a trajectoryvector, and 
														// returns the total length and trajectory vector (minus trajectory label).
};



#ifndef NO_INLINE
#include "inline.cpp"
#endif


#endif /*CORTICALSIM_H_*/
