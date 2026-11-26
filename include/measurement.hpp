#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include "types.hpp"
#include "parameters.hpp"

class Measurement
{
  public:

    double time;
    double length_density;
    double optical_density;
    double average_length;
    int number_of_m_ts;
    double segments_per_mt;
    OrderParameters order;
    int growing_number;
    int shrinking_number;
    int segments;
    int trajectories;
    int zipper_count;
    int crossover_count;
    int induced_catastrophe_count;
    int valid_d_event_count;
    int invalid_d_event_count;
    int s_event_count;
    int length_severing_count;
    int intersection_severing_count;
    int occupied_intersection_count;
    double g_eff_adjust;
};

ostream& operator<<(ostream&, const Measurement);
void write_measurement_descriptors(ostream&);

#endif // MEASUREMENT_HPP
