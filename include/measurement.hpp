#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include "types.hpp"
#include "parameters.hpp"

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

ostream& operator<<(ostream&, const Measurement);
void writeMeasurementDescriptors(ostream&);

#endif // MEASUREMENT_HPP
