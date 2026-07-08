#ifndef TYPES_HPP
#define TYPES_HPP

/*******************************************************************************************************************************************************/
// For details of simulation implementation:
// Chakrabortty B, Blilou I, Scheres B, Mulder BM (2018) A computational
// framework for cortical microtubule dynamics in realistically shaped plant
// cells. PLoS Comput Biol 14(2): e1005959.
// https://doi.org/10.1371/journal.pcbi.1005959

/***********************************************************************************************/
// This is a framework to simulate cortical microtubule dynamics in arbitrary
// Shaped plant cells. The program is developed by Bandan Chakrabortty
// (07/03/2013) and some basic implementations are used from (cortSim:1.20)
// developed by Simon Tindemans (01/04/2008)

/*********************************************************************************************/
/*
 * IMPORTANT: pitfalls to watch out for when extending the code
 * * The microtubule is only assumed to grow in length at the plus end. This is
 * not an MTTip property, but a property of the Microtubule functions.
 * Two-sided growth can be enabled, but needs plus/minus distinctions for the
 * events. Also, tips with velocity zero (usually non-treadmilling minus ends)
 * are assumed to be in a 'shrinking' state.
 * * The program assumes a stable sort order for the intersection lists. This
 * is true in practice (in all common implementations), but not mandated by the
 * standard.
 * * If the tubulin pool size is decreased on the fly, care should be taken
 * that the actual density doesn't exceed the pool size. In this situation, the
 * behaviour is undefined.
 * * In microtubule.cpp and system.cpp, compiler warnings are (should be)
 * issued for using the 'this' pointer within the initializer list, because the
 * 'this' pointer cannot be used until the initialization has completed. The
 *   current use is ok, because the pointer is only used to store the address
 * of the parent object. NOTE: when making changes to the code, take care NOT
 * to use this pointer within the constructors of other objects!
 */

#define PROGRAM_VERSION "BC.2018"
#define CROSS_SEV
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <boost/pool/pool_alloc.hpp>
#include <string>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "randhub.hpp"
#include "eig3.hpp"
#include "DLList.hpp"
#include <float.h>
#include "CompactList.hpp"
#include <sys/stat.h>
#include <cstddef>

using namespace Eigen;
using namespace std;

class System;
class Microtubule;
class Segment;
class MTTip;
class DeterministicEvent;
class Parameters;
class Geometry;
class Region;
class Trajectory;
class TrajectoryVector;
class Intersection;

#ifdef CROSS_SEV
class OccupiedIntersection;
#endif

typedef multimap<double, Intersection>::iterator IntersectionItr;
typedef int EventDescriptorIndex;
typedef int EventTrackingTag;
typedef list<Segment*>::iterator TrjSegmentTag;
typedef list<MTTip*>::iterator TrjMTTipTag;
typedef list<MTTip*>::iterator RegionMTTipTag;

// numeric parameters
const double PI = 3.141592653589793;
const double ZERO_CUTOFF = 1000000 * numeric_limits<double>::epsilon(); // approx 10E-10;//
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
typedef enum
{
    forward = 1,
    backward = -1
} Direction;

extern string DirectionTypeText[];

typedef enum
{
    ev_none,
    ev_wall = 1,
    ev_collision,
    ev_backtrack,
    ev_end_of_segment,
    ev_disappear
} DeterministicEventType;

typedef enum
{
    measure = 1,
    snapshot,
    stop,
    status,
    parameter_change,
    signalPPB
} GlobalEventType;

typedef enum
{
    catastrophe = 1,
    rescue,
    katanin,
    severingAtCross,
    nucleation
} StochasticEventType;

typedef enum
{
    ct_zipper = 1,
    ct_crossover,
    ct_inducedCatastrophe
} CollisionType;

typedef enum
{
    nuc_isotropic,
    nuc_COUNT_LAST
} NucleationType;

extern string NucleationTypeText[];

typedef enum
{
    int_zipFirst,
    int_catFirst,
    int_COUNT_LAST
} InteractionType;

extern string InteractionTypeText[];

typedef enum
{
    bdl_simple,
    bdl_sticky,
    bdl_noZip,
    bdl_multiCollision,
    bdl_Ncollision,
    bdl_COUNT_LAST
} BundleType;

extern string BundleTypeText[];

typedef enum
{
    t_minus,
    t_plus
} TipType;

typedef enum
{
    mt_growing,
    mt_shrinking
} MTType;

#endif // TYPES_HPP
