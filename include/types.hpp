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
enum class Direction
{
    FORWARD = 1,
    BACKWARD = -1
};

extern string direction_type_text[];

enum class DeterministicEventType : int
{
    NONE = 0,
    WALL = 1,
    COLLISION = 2,
    BACKTRACK = 3,
    END_OF_SEGMENT = 4,
    DISAPPEAR = 5
};

enum class GlobalEventType : int
{
    MEASURE = 1,
    SNAPSHOT = 2,
    STOP = 3,
    STATUS = 4,
    PARAMETER_CHANGE = 5,
    SIGNAL_PPB = 6
};

enum class StochasticEventType : int
{
    CATASTROPHE = 1,
    RESCUE = 2,
    KATANIN = 3,
    SEVERING_AT_CROSS = 4,
    NUCLEATION = 5
};

enum class CollisionType : int
{
    ZIPPER = 1,
    CROSSOVER = 2,
    INDUCED_CATASTROPHE = 3
};

enum class NucleationType : int
{
    ISOTROPIC = 0,
    COUNT_LAST = 1
};

extern string nucleation_type_text[];

enum class InteractionType : int
{
    ZIP_FIRST = 0,
    CAT_FIRST = 1,
    COUNT_LAST = 2
};

extern string interaction_type_text[];

enum class BundleType : int
{
    SIMPLE = 0,
    STICKY = 1,
    NO_ZIP = 2,
    MULTI_COLLISION = 3,
    NCOLLISION = 4,
    COUNT_LAST = 5
};

extern string bundle_type_text[];

enum class TipType : int
{
    MINUS = 0,
    PLUS = 1
};

enum class MTType : int
{
    GROWING = 0,
    SHRINKING = 1
};

#endif // TYPES_HPP
