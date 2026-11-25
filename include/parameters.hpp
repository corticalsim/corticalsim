#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <cstddef>
#include <numbers>
#include "types.hpp"
#include "linalg.hpp"

// numeric parameters
// const double PI = 3.141592653589793;
const double PI = std::numbers::pi_v<double>;
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
    map<string, double> edgCatMap;
    map<string, double> faceCatMap;
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
    void initialize(const char*);
    bool reinitialize(const char*);
    bool read_from_file(const char*, bool);
    bool write_to_file();
    void verify_parameters();
    bool calc_theory_parameters();
};

struct OrderParameters
{
    double R, C;
    double Rdirector[3];
    vector<double> localOrder;
    vector<Eigen::Vector3d> Sv;
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
    Eigen::Vector3d Sv;

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

        Sv << 0.0, 0.0, 0.0;

        return;
    }

    double extract_r(double director[3], string geometry)
    {
        Eigen::Matrix3d matrix{ { Qxx, Qxy, Qxz }, { Qxy, Qyy, Qyz }, { Qxz, Qyz, Qzz } };
        Eigen::Matrix3d evec_mat{ Eigen::Matrix3d::Zero() };
        Eigen::Vector3d e_val{ Eigen::Vector3d::Zero() };
        double sel_eigen_val{ 0.0 };

        eigen_decomposition(matrix, evec_mat, e_val);

        if (geometry == "2D-plane_1_0_0")
        {
            int max_pos(0);
            if (e_val[1] > e_val[0])
            {
                max_pos = 1;
            }
            if (e_val[2] > e_val[max_pos])
            {
                max_pos = 2;
            }
            sel_eigen_val = fabs(e_val[max_pos]);
            for (size_t i = 0; i < 3; ++i)
            {
                director[i] = evec_mat(i, max_pos);
            }
        }
        else
        {
            int min_pos(0);
            if (e_val[1] < e_val[0])
            {
                min_pos = 1;
            }
            if (e_val[2] < e_val[min_pos])
            {
                min_pos = 2;
            }
            sel_eigen_val = fabs(e_val[min_pos]);
            for (size_t i = 0; i < 3; ++i)
            {
                director[i] = evec_mat(i, min_pos);
            }
        }

        return (sel_eigen_val);
    }
};

#endif // PARAMETERS_HPP
