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
#include <float.h>
#include <sys/stat.h>
#include <cstddef>
#include "randhub.hpp"
#include "DLList.hpp"
#include "CompactList.hpp"

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

// new data type definition
typedef enum
{
    FORWARD = 1,
    BACKWARD = -1
} Direction;

extern string direction_type_text[];

typedef enum
{
    EV_NONE,
    EV_WALL = 1,
    EV_COLLISION,
    EV_BACKTRACK,
    EV_END_OF_SEGMENT,
    EV_DISAPPEAR
} DeterministicEventType;

typedef enum
{
    MEASURE = 1,
    SNAPSHOT,
    STOP,
    STATUS,
    PARAMETER_CHANGE,
    SIGNAL_PPB
} GlobalEventType;

typedef enum
{
    CATASTROPHE = 1,
    RESCUE,
    KATANIN,
    SEVERING_AT_CROSS,
    NUCLEATION
} StochasticEventType;

typedef enum
{
    CT_ZIPPER = 1,
    CT_CROSSOVER,
    CT_INDUCED_CATASTROPHE
} CollisionType;

typedef enum
{
    NUC_ISOTROPIC,
    NUC_COUNT_LAST
} NucleationType;

extern string nucleation_type_text[];

typedef enum
{
    INT_ZIP_FIRST,
    INT_CAT_FIRST,
    INT_COUNT_LAST
} InteractionType;

extern string interaction_type_text[];

typedef enum
{
    BDL_SIMPLE,
    BDL_STICKY,
    BDL_NO_ZIP,
    BDL_MULTI_COLLISION,
    BDL_NCOLLISION,
    BDL_COUNT_LAST
} BundleType;

extern string bundle_type_text[];

typedef enum
{
    T_MINUS,
    T_PLUS
} TipType;

typedef enum
{
    MT_GROWING,
    MT_SHRINKING
} MTType;

#endif // TYPES_HPP
