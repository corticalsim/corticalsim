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

    System* m_system;
    double m_v_plus;
    double m_v_min;
    double m_v_tm;
    double m_k_sev;
    double m_k_cross;
    double m_k_cat;
    double m_k_res;
    double m_k_nuc;
    double m_pp_bforming_time;
    double m_pool_density;
    int m_treadmilling_enabled;
    int m_severing_enabled;
    int m_cross_severing_enabled;
    int m_cross_severing_top;
    double m_cross_severing_start_angle;
    int m_ppb;
    int m_restricted_pool;
    int m_edge_catastrophe_enabled;
    int m_edge_catastrophe_smooth;
    double m_p_cat_special_edge_max;
    double m_p_cat_regular_edge_max;
    int m_edg_number;
    int m_face_number;
    map<string, double> m_edg_cat_map;
    map<string, double> m_face_cat_map;
    vector<double> m_region_kcat_multiplier;
    NucleationType m_nucleation_type;
    int m_zippering_enabled;
    int m_catastrophes_enabled;
    int m_proportional_catastrophes;
    double m_induced_catastrophe_fraction;
    double m_zip_fraction;
    double m_pp_bk_nuc_fraction;
    double m_cat_start_angle;
    double m_magic_angle;
    InteractionType m_interaction_type;
    BundleType m_bundle_type;
    string m_geom_param;
    unsigned long m_seed;
    double m_stop_time;
    double m_measurement_interval;
    double m_wall_clock_limit;
    double m_memory_limit;
    string m_input_dir;
    string m_output_dir;
    string m_movie_dir;
    int m_create_subdir;
    double m_new_parameter_read_interval;
    string m_new_parameter_file;
    int m_movie_enabled;
    int m_show_output;
    int m_show_mesh;
    string m_geometry;
    double m_movie_frame_interval;
    double m_geometry_scale_factor;

    double m_c0calc;
    double m_x0calc;
    double m_z0calc;

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
    double m_r, m_c;
    double m_rdirector[3];
    vector<double> m_local_order;
    vector<Eigen::Vector3d> m_sv;
};

class OrderParametersRaw
{
  public:

    double m_local_l;
    double m_local_l_opt;
    double m_qxx;
    double m_qxy;
    double m_qxz;
    double m_qyy;
    double m_qyz;
    double m_qzz;
    double m_local_order;
    Eigen::Vector3d m_sv;

    OrderParametersRaw(void)
    {
        m_local_l = 0.;
        m_local_l_opt = 0.;
        m_qxx = 0;
        m_qxy = 0;
        m_qxz = 0;
        m_qyy = 0;
        m_qyz = 0;
        m_qzz = 0;

        m_sv << 0.0, 0.0, 0.0;

        return;
    }

    double extract_r(double director[3], string geometry)
    {
        Eigen::Matrix3d matrix{ { m_qxx, m_qxy, m_qxz }, { m_qxy, m_qyy, m_qyz }, { m_qxz, m_qyz, m_qzz } };
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
