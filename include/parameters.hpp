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
    double v_plus;
    double v_min;
    double v_tm;
    double k_sev;
    double k_cross;
    double k_cat;
    double k_res;
    double k_nuc;
    double pp_bforming_time;
    double pool_density;
    int treadmilling_enabled;
    int severing_enabled;
    int cross_severing_enabled;
    int cross_severing_top;
    double cross_severing_start_angle;
    int ppb;
    int restricted_pool;
    int edge_catastrophe_enabled;
    int edge_catastrophe_smooth;
    double p_cat_special_edge_max;
    double p_cat_regular_edge_max;
    int edg_number;
    int face_number;
    map<string, double> edg_cat_map;
    map<string, double> face_cat_map;
    vector<double> region_kcat_multiplier;
    NucleationType nucleation_type;
    int zippering_enabled;
    int catastrophes_enabled;
    int proportional_catastrophes;
    double induced_catastrophe_fraction;
    double zip_fraction;
    double pp_bk_nuc_fraction;
    double cat_start_angle;
    double magic_angle;
    InteractionType interaction_type;
    BundleType bundle_type;
    string geom_param;
    unsigned long seed;
    double stop_time;
    double measurement_interval;
    double wall_clock_limit;
    double memory_limit;
    string input_dir;
    string output_dir;
    string movie_dir;
    int create_subdir;
    double new_parameter_read_interval;
    string new_parameter_file;
    int movie_enabled;
    int show_output;
    int show_mesh;
    string geometry;
    double movie_frame_interval;
    double geometry_scale_factor;

    double c0_calc;
    double x0_calc;
    double z0_calc;

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
    double r, c;
    double rdirector[3];
    vector<double> local_order;
    vector<Eigen::Vector3d> sv;
};

class OrderParametersRaw
{
  public:

    double local_l;
    double local_l_opt;
    double qxx;
    double qxy;
    double qxz;
    double qyy;
    double qyz;
    double qzz;
    double local_order;
    Eigen::Vector3d sv;

    OrderParametersRaw(void)
    {
        local_l = 0.;
        local_l_opt = 0.;
        qxx = 0;
        qxy = 0;
        qxz = 0;
        qyy = 0;
        qyz = 0;
        qzz = 0;

        sv << 0.0, 0.0, 0.0;

        return;
    }

    double extract_r(double director[3], string geometry)
    {
        Eigen::Matrix3d matrix{ { qxx, qxy, qxz }, { qxy, qyy, qyz }, { qxz, qyz, qzz } };
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
