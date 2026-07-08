#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "types.hpp"

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
    bool readFromFile(const char*, bool);
    bool writeToFile();
    void verifyParameters();
    bool calcTheoryParameters();
};

struct OrderParameters
{
    double R, C;
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

        Sv << 0.0, 0.0, 0.0;

        return;
    }

    double extractR(double director[3], string geometry)
    {
        double matrix[3][3] = { { Qxx, Qxy, Qxz }, { Qxy, Qyy, Qyz }, { Qxz, Qyz, Qzz } };
        double evecMat[3][3];
        double selEigenVal, eVal[3];
        eigen_decomposition(matrix, evecMat, eVal);

        if (geometry == "2D-plane_1_0_0")
        {
            int maxPos(0);
            if (eVal[1] > eVal[0])
            {
                maxPos = 1;
            }
            if (eVal[2] > eVal[maxPos])
            {
                maxPos = 2;
            }
            selEigenVal = fabs(eVal[maxPos]);
            for (int i = 0; i < 3; i++)
            {
                director[i] = evecMat[i][maxPos];
            }
        }
        else
        {
            int minPos(0);
            if (eVal[1] < eVal[0])
            {
                minPos = 1;
            }
            if (eVal[2] < eVal[minPos])
            {
                minPos = 2;
            }
            selEigenVal = fabs(eVal[minPos]);
            for (int i = 0; i < 3; i++)
            {
                director[i] = evecMat[i][minPos];
            }
        }

        return (selEigenVal);
    }
};

#endif // PARAMETERS_HPP
