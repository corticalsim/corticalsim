#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include "types.hpp"
#include "parameters.hpp"

class Measurement
{
  public:

    double m_time;
    double m_length_density;
    double m_optical_density;
    double m_average_length;
    int m_number_of_m_ts;
    double m_segments_per_mt;
    OrderParameters m_order;
    int m_growing_number;
    int m_shrinking_number;
    int m_segments;
    int m_trajectories;
    int m_zipper_count;
    int m_crossover_count;
    int m_induced_catastrophe_count;
    int m_valid_d_event_count;
    int m_invalid_d_event_count;
    int m_s_event_count;
    int m_length_severing_count;
    int m_intersection_severing_count;
    int m_occupied_intersection_count;
    double m_g_eff_adjust;
};

ostream& operator<<(ostream&, const Measurement);
void write_measurement_descriptors(ostream&);

#endif // MEASUREMENT_HPP
