#include "mesh.hpp"
#include "region.hpp"
#include "system.hpp"
#include "mt_tip.hpp"
#include "parameters.hpp"

void System::collision_probabilities(double angle, double& pcat, double& pzip)
{
    switch (m_p.m_interaction_type)
    {
        case INT_ZIP_FIRST:

            // angle < magic angle
            if (angle < PI * m_p.m_magic_angle / 180.0)
            {
                pcat = 0.0;

                // zippering enabled
                if (m_p.m_zippering_enabled)
                {
                    pzip = m_p.m_zip_fraction;
                }
                else
                {
                    pzip = 0.0;
                }
                break;
            }

            pzip = 0.0;

            // angle > p.catStartAngle and induced catastrophes enabled
            if (m_p.m_catastrophes_enabled && (angle > m_p.m_cat_start_angle / 180.0 * PI))
            {
                // proportional catastrophe enabled
                if (m_p.m_proportional_catastrophes)
                {
                    pcat = (angle - m_p.m_cat_start_angle / 180.0 * PI) / (PI / 2. - m_p.m_cat_start_angle / 180.0 * PI)
                           * m_p.m_induced_catastrophe_fraction;
                }

                // proportional catastrophe disabled
                else
                {
                    pcat = m_p.m_induced_catastrophe_fraction;
                }
            }
            else
            {
                pcat = 0.0;
            }

            break;

        case INT_CAT_FIRST:

            // angle > p.catStartAngle and induced catastrophes enabled
            if (m_p.m_catastrophes_enabled && (angle > m_p.m_cat_start_angle / 180.0 * PI))
            {
                // proportional catastrophe enabled
                if (m_p.m_proportional_catastrophes)
                {
                    pcat = (angle - m_p.m_cat_start_angle / 180.0 * PI) / (PI / 2. - m_p.m_cat_start_angle / 180.0 * PI)
                           * m_p.m_induced_catastrophe_fraction;
                }

                // proportional catastrophe disabled
                else
                {
                    pcat = m_p.m_induced_catastrophe_fraction;
                }
            }
            else
            {
                pcat = 0.0;
            }

            // angle < magic angle and zippering enabled
            if (m_p.m_zippering_enabled && (angle < m_p.m_magic_angle / 180.0 * PI))
            {
                pzip = 1.0 - pcat;
            }
            else
            {
                pzip = 0.0;
            }
            break;

        default:
            break;
    }

#ifdef DBG_INTERACTION_TEST
    cout << "InteractionTest\t" << angle << "\t" << Pzip << "\t" << Pcat << "\t" << 1 - Pzip - Pcat << "\n";
#endif

    return;
}

double System::multi_pcross(int npar, int ncoll, double x_single, double z_single)
{
    double total_sum = 0;

    int i, j;

    for (i = 0; i <= npar; i++)
    {
        for (j = max(0, i - ncoll); j <= max(0, i - 1); j++)
        {
            total_sum += (npar + 1 - i) * m_binomial_table[ncoll][i - j] * m_binomial_table[max(0, i - 1)][j]
                        * pow(x_single, 2 * j + ncoll - i) * pow(z_single, 2 * (i - j));
        }
    }
    return total_sum / (npar + 1);
}

double System::multi_pzip(int n_par, int n_coll, double x_single, double z_single)
{
    double total_sum = 0;

    int i, j, k;

    for (i = 0; i <= n_par; i++)
    {
        for (j = 0; j <= n_coll - 1; j++)
        {
            for (k = max(0, i - j); k <= i; k++)
            {
                total_sum += m_binomial_table[j][i - k] * m_binomial_table[i][k] * pow(x_single, 2 * k + j - i)
                            * pow(z_single, 2 * (i - k) + 1);
            }
        }
    }
    return total_sum / (n_par + 1);
}

CollisionType System::collision_result(double angle, int n_parallel, int n_cross) // range: 0..PI/2
{
#ifdef DBG_ASSERT
    if ((angle < 0) || (angle > PI / 2))
    {
        cout << "ERROR: angle out of range [0..PI/2] : " << angle << "\n";
        exit(-2);
    }
#endif

#ifdef DBG_ASSERT
    {
        if (nCross == 0)
        {
            cerr << "ERROR: crossover collision passed to collisionResult()\n";
            return ct_crossover;
        }
    }
#endif

    double pcat{ 0.0 };
    double pzip{ 0.0 };
    double psum{ 0.0 };

    // get a random number
    double rand = m_random_gen.rand_dbl_exc();

    // get the available Pcat and Pzip
    collision_probabilities(angle, pcat, pzip);

    switch (m_p.m_bundle_type)
    {
        // TODO
        // This needs to be fixed properly.
        // Adding this case or `default` at the end causes an error becasue
        // of the initialisation of multiPx & multiPz in the last case.
        case BDL_COUNT_LAST:
            break;
        case BDL_NO_ZIP:

            if (n_parallel > 0)
            {
                if (rand < pcat)
                {
                    return CT_INDUCED_CATASTROPHE;
                }
            }
            // break omitted on purpose!

        case BDL_STICKY:

            if (n_parallel > 0)
            {
                return CT_CROSSOVER;
            }
            // break omitted on purpose!

        case BDL_SIMPLE:

            // induced catastrophe selected
            if (rand < pcat)
            {
                return CT_INDUCED_CATASTROPHE;
            }

            // zippering selected
            else if (rand < pcat + pzip)
            {
                return CT_ZIPPER;
            }

            // cross over slected
            else
            {
                return CT_CROSSOVER;
            }

            // though never reached
            // break;

        case BDL_NCOLLISION:

            // add up total probabilities
            pcat *= n_cross;
            pzip *= n_cross;
            psum = pcat + pzip;

            // now normalize the probabilities
            if (psum > 1)
            {
                pcat /= psum;
                pzip /= pzip;
            }

            if (rand < pcat)
            {
                return CT_INDUCED_CATASTROPHE;
            }

            else if (rand < pcat + pzip)
            {
                return CT_ZIPPER;
            }

            else
            {
                return CT_CROSSOVER;
            }

            // though never reached
            // break;

        case BDL_MULTI_COLLISION:

            double pcross = 1.0 - pcat - pzip;
            double multi_px = multi_pcross(min(n_parallel, MAXBINOM - 1), min(n_cross, MAXBINOM - 1), pcross, pzip);
            double multi_pz = multi_pzip(min(n_parallel, MAXBINOM - 1), min(n_cross, MAXBINOM - 1), pcross, pzip);

            if (rand < multi_px)
            {
                return CT_CROSSOVER;
            }

            else if (rand < multi_px + multi_pz)
            {
                return CT_ZIPPER;
            }

            else
            {
                return CT_INDUCED_CATASTROPHE;
            }

            // though never reached
            // break;
    }

    // should never occur
    return CT_CROSSOVER;
}

System::System(char* par_file):
    m_geometry(nullptr),
    m_p(this),
    m_system_time(0),
    m_system_time_offset(0),
    m_current_time_tag(0),
    m_total_length(0),
    m_mem_usage(0),
    m_stop_signal(false),
    m_wall_clock_start_time(time(nullptr)),
    m_total_s_event_count(0),
    m_total_valid_d_event_count(0),
    m_total_invalid_d_event_count(0),
    m_total_length_severing_count(0),
    m_total_intersection_severing_count(0),
    m_boundary_crossing_count(0),
    m_total_zipper_count(0),
    m_total_crossover_count(0),
    m_total_induced_catastrophe_count(0),
    m_count_segments(0),
    m_count_trajectories(0),
    m_count_intersections(0),
    m_time_queue(this, &System::identity, &System::identity),
    m_v_plus_queue(this, &System::v_plus_to_time, &System::time_to_v_plus),
    m_next_status_event_time(0),
    m_next_snapshot_event_time(0),
    m_event_descriptor_id(0),
    m_event_id(0),
    m_data_saved(false)
{
    /// read parameters
    m_p.initialize(par_file);

    // initialize number of growing tips in different faces
    for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
    {
        m_growing_tips_reg.push_back(0);
    }

    // load geometry
    m_geometry = new TriMeshGeometry(this);

    /// initialize the random seed with the desired one
    m_random_gen.seed(m_p.m_seed);

    /// used in multiPcross & multiPzip events
    make_binomial_table();

    /// perform further initialization
    initialize_output();
    m_p.write_to_file();

    /// start creating MT events
    determine_stochastic_event();
    m_time_queue.push_global(m_p.m_measurement_interval, STATUS);
    m_next_status_event_time = m_p.m_measurement_interval;

    /// if snapShot events are requested
    if (m_p.m_movie_enabled)
    {
        m_time_queue.push_global(m_p.m_movie_frame_interval, SNAPSHOT);
        m_next_snapshot_event_time = m_p.m_movie_frame_interval;
    }

    /// if parameter file is updated on run
    if (m_p.m_ppb)
    {
        m_time_queue.push_global(m_p.m_pp_bforming_time, SIGNAL_PPB);
        m_next_ppb_event_time = m_p.m_pp_bforming_time;
    }

    /// if parameter file is updated on run
    if (m_p.m_new_parameter_read_interval > ZERO_CUTOFF)
    {
        m_time_queue.push_global(m_p.m_new_parameter_read_interval, PARAMETER_CHANGE);
        m_next_parameter_event_time = m_p.m_new_parameter_read_interval;
    }

    return;
}

System::~System()
{
    // destroy all MTs
    m_growing_mts.remove_all();
    m_shrinking_mts.remove_all();

    return;
}

void System::run(double sim_time, string& c_s, string& m_s)
{
    int wall_time_polling_counter(0);
    int memory_polling_counter(0);
    m_time_queue.push_global(sim_time - m_system_time_offset, STOP);

    if (m_p.m_show_output == 1)
    {
        cout << "===============================================" << endl;
        cout << " Time stamped simulation output " << endl;
        cout << "===============================================" << endl;
    }

    while (!m_stop_signal)
    {
        if (++wall_time_polling_counter == CLOCK_POLLING_INTERVAL)
        {
            if (difftime(time(nullptr), m_wall_clock_start_time) >= m_p.m_wall_clock_limit)
            {
                // time's up; send a stop signal
                cerr << "Wall clock limit (" << m_p.m_wall_clock_limit << " seconds) reached. Stopping.\n";
                m_time_queue.push_global(m_system_time, STOP);

                // write final snapshot when running out of time.
                if (m_p.m_movie_enabled)
                {
                    m_geometry->output_snapshot(m_movie_file);
                }

                // emergency exit the program
                emergency_break();
                exit(-1);
            }
            wall_time_polling_counter = 0;
        }

        if (++memory_polling_counter == MEMORY_POLLING_INTERVAL)
        {
            if (estimate_memory_footprint() > 1024 * 1024 * m_p.m_memory_limit)
            {
                cerr << "Memory limit (approx. " << m_p.m_memory_limit << " MB) reached. Stopping.\n";
                m_time_queue.push_global(m_system_time, STOP);

                // emergency exit the program
                emergency_break();
                exit(-1);
            }
            memory_polling_counter = 0;
        }

        if (m_event_id > QUEUE_FLUSH_INTERVAL)
        {
#ifdef DBG_QMGMT
            cout << "Flushing queues to reduce numerical drift.\n";
#endif

            // moved forward (used to be after flush and reload)
            if (!integrity_check())
            {
                exit(-1);
            }

            flush_and_reload();
            if (!integrity_check())
            {
                exit(-1);
            }

#ifdef DBG_QMGMT
            cout << "Queue reloaded.\n";
#endif
        }

        // continue to next event
        next_event();
    }

    // stop signal received : write down necessary output(s)
    write_measurements_to_file(0);
    close_files();

    return;
}

void System::emergency_break()
{
    // write current output(s) and close all opened files
    write_measurements_to_file(0);
    close_files();

    return;
}

int System::estimate_memory_footprint()
{
    // calculate memory in use
    int total(0);
    total += (m_growing_mts.size() + m_shrinking_mts.size()) * (sizeof(Microtubule) + 2 * sizeof(MTTip*));
    total += m_event_descriptor_map.size() * (sizeof(pair<EventDescriptorIndex, EventDescriptor*>) + 3 * sizeof(void*));
    total += m_count_segments * (sizeof(Segment) + sizeof(Segment*));
    total += m_count_trajectories * (sizeof(Trajectory));
    total += m_count_intersections * (2 * (sizeof(pair<double, Intersection>) + 3 * sizeof(void*)));
    total += (m_time_queue.m_queue.size() + m_v_plus_queue.m_queue.size()) * (sizeof(DeterministicEvent) + 2 * sizeof(void*));
    total += m_measurement_history.size() * sizeof(Measurement);

    return total;
}

bool System::integrity_check()
{
    bool valid(true);

    // integrity check of the gemoetry embeded components
    if (!m_geometry->integrity_check())
    {
        valid = false;
    }

    // integrity check of each growing MTs
    Microtubule* mt = m_growing_mts.first();
    while (mt != nullptr)
    {
        if (!mt->integrity_check())
        {
            valid = false;
        }

        mt = mt->next();
    }

    // integrity check of each shrinking MTs
    mt = m_shrinking_mts.first();
    while (mt != nullptr)
    {
        if (!mt->integrity_check())
        {
            valid = false;
        }

        mt = mt->next();
    }

    return valid;
}

void System::flush_and_reload(bool refresh_parameter_event)
{
    double time_offset(m_system_time);

    // update all microtubule lengths
    update_all();

    // recalculate total lengths to avoid drift (in DBG mode, also check for large deviations)
    double region_length(0.);
    double system_length(0);
    Trajectory* tptr(nullptr);

    for (size_t ridx = 0; ridx < m_geometry->m_regions.size(); ridx++)
    {
        region_length = 0;
        tptr = m_geometry->m_regions[ridx]->m_trajectories.first();

        while (tptr != nullptr)
        {
            region_length += tptr->segment_length();
            tptr = tptr->next();
        }

        if (abs(m_geometry->m_regions[ridx]->m_total_length - region_length) > ZERO_CUTOFF * max(100., region_length))
        {
            cerr << "Unacceptable drift in region length (drift = "
                 << abs(m_geometry->m_regions[ridx]->m_total_length - region_length)
                 << "). Exiting. [try lowering flush interval]\n";
            emergency_break();
            exit(-1);
        }
        m_geometry->m_regions[ridx]->m_total_length = region_length;
        system_length += region_length;
    }

    if (abs(m_total_length - system_length) > ZERO_CUTOFF * max(1., m_total_length))
    {
        cerr << "Unacceptable drift in system length (drift = " << abs(m_total_length - system_length)
             << "). Exiting. [try lowering flush interval]\n";
        emergency_break();
        exit(-1);
    }

    m_total_length = system_length;

    // reset master timer
    m_system_time_offset += time_offset;

    m_system_time = 0;
    m_time_queue.flush();
    m_v_plus_queue.flush();
    m_event_descriptor_map.clear();
    m_event_id = 0;
    m_event_descriptor_id = 0;

    // now, reload all events (this is a PROBLEM for global events that have not been saved. They need to be re-inserted
    // manually)
    determine_stochastic_event();
    m_time_queue.push_global(m_p.m_stop_time - m_system_time_offset, STOP);
    m_next_status_event_time -= time_offset;
    m_time_queue.push_global(m_next_status_event_time, STATUS);

    if (refresh_parameter_event)
    {
        m_next_parameter_event_time -= time_offset;
        m_time_queue.push_global(m_next_parameter_event_time, PARAMETER_CHANGE);
    }

    if (m_p.m_movie_enabled)
    {
        m_next_snapshot_event_time -= time_offset;
        m_time_queue.push_global(m_next_snapshot_event_time, SNAPSHOT);
    }

    if (m_p.m_ppb)
    {
        m_next_ppb_event_time -= time_offset;
        m_time_queue.push_global(m_next_ppb_event_time, SIGNAL_PPB);
    }

    Microtubule* mt(nullptr);
    Segment* seg(nullptr);
    mt = m_growing_mts.first();

    while (mt != nullptr)
    {
        // this can go into a change time function
        mt->m_nucleation_time -= time_offset;
        seg = mt->m_segments.first();

        while (seg != nullptr)
        {
            seg->m_nucleation_time -= time_offset;
            seg = seg->next();
        }

        mt->m_plus.m_event.m_index = register_event_descriptor(&(mt->m_plus.m_event));
        mt->m_minus.m_event.m_index = register_event_descriptor(&(mt->m_minus.m_event));
        mt->m_disappear_event.m_index = register_event_descriptor(&(mt->m_disappear_event));
        mt->m_plus.determine_event();
        mt->m_minus.determine_event();
        mt->set_disappear_event();
        mt = mt->next();
    }

    mt = m_shrinking_mts.first();

    while (mt != nullptr)
    {
        mt->m_nucleation_time -= time_offset;
        seg = mt->m_segments.first();

        {
            seg->m_nucleation_time -= time_offset;
            seg = seg->next();
        }

        mt->m_plus.m_event.m_index = register_event_descriptor(&(mt->m_plus.m_event));
        mt->m_minus.m_event.m_index = register_event_descriptor(&(mt->m_minus.m_event));
        mt->m_disappear_event.m_index = register_event_descriptor(&(mt->m_disappear_event));
        mt->m_plus.determine_event();
        mt->m_minus.determine_event();
        mt->set_disappear_event();
        mt = mt->next();
    }

    return;
}

double lambert_w(double x)
{
    double out(0.0);
    const double exp_min1 = exp(-1.0);
    const double zero_cutoff = pow(2.0, -DBL_MANT_DIG / 6.0);
    const double eta_cutoff = 34.0 * pow(2.0, 2.0 * (1.0 - DBL_MANT_DIG) / 7.0);

    // LambertW function (Barry et al., ACM Trans. on Math. Software, Vol. 21, No. 2, June 1995, pp. 161--171)
    if (x < -exp_min1)
    {
        return -VERY_LARGE;
    }

    if (abs(x) < zero_cutoff)
    {
        out = x * (60 + 114 * x + 17 * x * x) / (60 + 174 * x + 101 * x * x);
    }
    else if (x > 20.0)
    {
        double h = exp(-1.124491989777808 / (0.4225028202459761 + log(x)));
        out = log(x / log(x / pow(log(x), h)));

        double z_n = log(x / out) - out;
        double temp = 2.0 * (1.0 + out) * (1.0 + out + z_n * 2.0 / 3.0);
        double e_n = (z_n / (1.0 + out)) * (temp - z_n) / (temp - 2.0 * z_n);
        out = out * (1.0 + e_n);
    }
    else
    {
        double eta = 2.0 + 2.0 * exp(1.0) * x;
        double sqrt_eta = sqrt(eta);
        double d(0.0);
        if (sqrt_eta >= eta_cutoff)
        {
            double n2 = 4.612634277343749 * sqrt(sqrt(sqrt_eta + 1.09556884765625));
            double n1 = (4.0 - 3.0 * sqrt(2.0) + (2.0 * sqrt(2.0) - 3.0) * n2) / (sqrt(2.0) - 2.0);
            d = n1 * sqrt_eta / (n2 + sqrt_eta);

            out = -1.0 + sqrt_eta / (1.0 + sqrt_eta / (3.0 + d));

            double z_n = log(x / out) - out;
            double temp = 2.0 * (1.0 + out) * (1.0 + out + z_n * 2.0 / 3.0);
            double e_n = (z_n / (1.0 + out)) * (temp - z_n) / (temp - 2.0 * z_n);
            out = out * (1.0 + e_n);
        }
        else
        {
            d = sqrt_eta
                / (8.0 / 3.0 + sqrt_eta / (135.0 / 83.0 + sqrt_eta / (166.0 / 39.0 + sqrt_eta * 3167.0 / 3549.0)));
            out = -1.0 + sqrt_eta / (1.0 + sqrt_eta / (3.0 + d));
        }
    }

#ifdef DBG_ACID_TEST
    double test = (x - out * exp(out)) / ((1.0 + out) * exp(out));
    if (abs(test) > ZERO_CUTOFF * abs(x))
    {
        cerr << "LambertW() accuracy [" << test << "] below threshold. Output may contain inaccuracies.\n";
    }
#endif

    return out;
}

double System::v_plus_to_time(double dist)
{
    // for infinite pool (implimented as identity function)
    if (m_p.m_restricted_pool == 0)
    {
        return dist;
    }

    // for finite tubulin pool
    double lmax = (m_p.m_pool_density * m_geometry->m_area);
    double inv_l_max = 1.0 / (m_p.m_pool_density * m_geometry->m_area);

    double alpha = (m_p.m_v_plus - m_p.m_v_tm) * m_growing_mts.size() + (m_p.m_v_min - m_p.m_v_tm) * m_shrinking_mts.size();

    if (m_growing_mts.size() == 0)
    {
        if (m_shrinking_mts.size() == 0)
        {
            return dist / (1.0 - m_total_length * inv_l_max);
        }
        else
        {
            return ((lmax - m_total_length) - sqrt(pow(lmax - m_total_length, 2) - 2. * alpha * lmax * dist)) / alpha;
        }
    }

    double beta = m_p.m_v_plus * m_growing_mts.size() * inv_l_max;

    if ((m_shrinking_mts.size() == 0) && ((m_p.m_v_tm / m_p.m_v_plus) < ZERO_CUTOFF))
    {
        return -log(1.0 - m_growing_mts.size() * m_p.m_v_plus * dist / (lmax - m_total_length)) / beta;
    }

    double pv = (alpha / beta - m_total_length) / (lmax - alpha / beta);
    double qv = m_growing_mts.size() * m_p.m_v_plus * dist / (lmax - alpha / beta);

    // verify whether this is an optimal trade-off
    if (pv - qv < 250)
    {
        return (qv - pv + lambert_w(pv * exp(pv - qv))) / beta;
    }

    // if not, do asymptotic approximation
    double l1 = pv - qv + log(pv);
    double l2 = log(l1);
    return (1.0 / beta)
           * (qv - pv + l1 - l2 + l2 / l1 + l2 * (l2 - 2) / (2 * pow(l1, 2))
              + l2 * (6 - 9 * l2 + 2 * pow(l2, 2)) / (6 * pow(l1, 3))
              + l2 * (-12 + 36 * l2 - 22 * pow(l2, 2) + 3 * pow(l2, 3)) / (12 * pow(l1, 4)));
}

double System::time_to_v_plus(double time)
{
    // for infinite pool (implimented as identity function)
    if (m_p.m_restricted_pool == 0)
    {
        return time;
    }

    // for finite tubulin pool
    double inv_l_max = 1.0 / (m_p.m_pool_density * m_geometry->m_area);
    double alpha = (m_p.m_v_plus - m_p.m_v_tm) * m_growing_mts.size() + (m_p.m_v_min - m_p.m_v_tm) * m_shrinking_mts.size();

    if (m_growing_mts.size() == 0)
    {
        return ((1.0 - m_total_length * inv_l_max) - 0.5 * alpha * inv_l_max * time) * time;
    }

    double beta = m_p.m_v_plus * m_growing_mts.size() * inv_l_max;
    double lambda_inf = inv_l_max * alpha / beta;
    double lambda_0 = m_total_length * inv_l_max;

    return (1.0 - lambda_inf) * time - (lambda_0 - lambda_inf) * (1.0 - exp(-beta * time)) / beta;
}

void System::advance_time(double new_time)
{
#ifdef DBG_SYSTEM
    cout << "DBG/SYSTEM: advancing time by " << newTime - systemTime << " seconds. New system time is " << newTime
         << ".\n";
#endif

    // make a jump in time step (taking advantage of using event-driven algorithm)
    m_time_queue.advance_time(new_time);
    m_v_plus_queue.advance_time(new_time);
    m_system_time = new_time;
    m_total_length += m_time_queue.progression(m_current_time_tag)
                   * (m_p.m_v_min * m_shrinking_mts.size() - m_p.m_v_tm * (m_growing_mts.size() + m_shrinking_mts.size()))
                   + m_v_plus_queue.progression(m_current_time_tag) * (m_p.m_v_plus * m_growing_mts.size());

    if (++m_current_time_tag == POSITION_CACHE_SIZE)
    {
        m_current_time_tag = 0;
        update_all(true);
    }

    m_time_queue.store_time(m_current_time_tag);
    m_v_plus_queue.store_time(m_current_time_tag);

    return;
}

void System::update_all(bool force_update)
{
    // always nice to be updated in lengths of all regions and microtubules (convenient for snapshots, etc.)
    for (size_t ridx = 0; ridx < m_geometry->m_regions.size(); ridx++)
    {
        m_geometry->m_regions[ridx]->update_region_length(force_update);
    }

    Microtubule* mt;
    mt = m_growing_mts.first();

    while (mt != nullptr)
    {
        mt->update_length(force_update);

#ifdef DBG_ACID_TEST
        if (mt->segments.size() == 1)
        {
            if ((mt->segments.first()->end - mt->segments.first()->start) * mt->segments.first()->dir < -ZERO_CUTOFF)
            {
                cerr << "ERROR: negative segment length after updating.\n";
                exit(-1);
            }
        }
#endif

        mt = mt->next();
    }

    mt = m_shrinking_mts.first();
    while (mt != nullptr)
    {
        mt->update_length(force_update);

#ifdef DBG_ACID_TEST
        if (mt->segments.size() == 1)
        {
            if ((mt->segments.first()->end - mt->segments.first()->start) * mt->segments.first()->dir < -ZERO_CUTOFF)
            {
                cerr << "ERROR: negative segment length after updating.\n";
                exit(-1);
            }
        }
#endif
        mt = mt->next();
    }

    return;
}

void System::next_event()
{
    double tq_value = numeric_limits<double>::max();
    double dq_value = numeric_limits<double>::max();

    // maybe clean up a little
    if (!m_time_queue.empty())
    {
        tq_value = m_time_queue.first_event_time();
    }

    if (!m_v_plus_queue.empty())
    {
        dq_value = m_v_plus_queue.first_event_time();
    }

    if ((m_next_stochastic_event_time < tq_value) && (m_next_stochastic_event_time < dq_value))
    {
        // execute stochastic event
        m_total_s_event_count++;
#ifdef DBG_EVENT
        cout << "DBG/Event: executing stochastic event [time=" << nextStochasticEventTime
             << ", type=" << nextStochasticEventType << "]\n";
#endif

        advance_time(m_next_stochastic_event_time);

        // execute next stochastic event
        switch (m_next_stochastic_event_type)
        {
            case CATASTROPHE:
                handle_catastrophe_event();
                break;
            case RESCUE:
                handle_rescue_event();
                break;
            case KATANIN:
                handle_severing_event();
                break;
#ifdef CROSS_SEV
            case SEVERING_AT_CROSS:
                handle_severing_at_cross_event();
                break;
#endif
            case NUCLEATION:
                handle_nucleation_event();
                break;
        }
        determine_stochastic_event();
    }

    else
    {
        double next_event_time(0.0);

        // execute deterministic event
        DeterministicEvent event;
        if (tq_value < dq_value)
        {
            event = m_time_queue.pop();
            next_event_time = tq_value;
        }
        else
        {
            event = m_v_plus_queue.pop();
            next_event_time = dq_value;
        }

        if (event.m_info_idx != -1)
        {
            map<EventDescriptorIndex, EventDescriptor*>::iterator event_itr = m_event_descriptor_map.find(event.m_info_idx);

            // whether the tip still exists and the tag still matches (valid event)
            if ((event_itr != m_event_descriptor_map.end()) && (event_itr->second->m_tag == event.m_tag))
            {
#ifdef DBG_EVENT
                cout << "DBG/Event: executing deterministic event [time=" << event.eventTimeDist
                     << ", type=" << eventItr->second->type << ", tag=" << event.tag << "]\n";
#endif

                m_total_valid_d_event_count++;
                advance_time(next_event_time);
                event_itr->second->m_mt->handle_event(event_itr->second);
                determine_stochastic_event();
            }

            // invalid event
            else
            {
                m_total_invalid_d_event_count++;
#ifdef DBG_EVENT
                cout << "DBG/Event: Event is no longer valid. [tag= " << event.tag << "]\n";
#endif
            }
        }

        // execute global event
        else
        {
#ifdef DBG_EVENT
            cout << "DBG/Event: executing global event [time=" << event.eventTimeDist << ", type=" << event.global_type
                 << "]\n";
#endif

            advance_time(event.m_event_time_dist);
            handle_global_event(event);
        }
    }

    return;
}

void System::handle_catastrophe_event()
{

    size_t tip_number(0), total_relevant_tips(0), ridx(0);
    int domain_type(0);
    RegionMTTipTag tip_tag;

    double norm(0.0);
    for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
    {
        norm += m_growing_tips_reg[d_tag] * m_p.m_region_kcat_multiplier[d_tag];
    }

    double rand_gen = m_random_gen.rand();

    double norm_fraction(0.0);
    for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
    {
        norm_fraction += m_growing_tips_reg[d_tag] * m_p.m_region_kcat_multiplier[d_tag];
        if (rand_gen < (double)(norm_fraction) / norm)
        {
            domain_type = d_tag;
            total_relevant_tips = m_growing_tips_reg[d_tag];
            break;
        }
    }

    tip_number = m_random_gen.rand_int(total_relevant_tips - 1);

    // determine the region this tip is located in
    // depending on whether it's a high or low number, start scanning from the beginning or end
    if (tip_number < total_relevant_tips / 2)
    {
        while (true)
        {
            if (m_geometry->m_regions[ridx]->m_face_tag == domain_type)
            {
                if (tip_number < m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size())
                {
                    break;
                }

                tip_number -= m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size();
            }
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tip_number = total_relevant_tips - 1 - tip_number;
        ridx = m_geometry->m_regions.size() - 1;
        while (true)
        {
            if (m_geometry->m_regions[ridx]->m_face_tag == domain_type)
            {
                if (tip_number < m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size())
                {
                    break;
                }
                tip_number -= m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size();
            }
            ridx--;
        }
        // flip the direction of the result for the next step
        tip_number = m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size() - 1 - tip_number;
    }

    // Now, select the tip with # 'tipNumber' from the current region.
    // Again, if the number is high within the number of tips in the region, start from the end
    if (tip_number < m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size() / 2)
    {
        tip_tag = m_geometry->m_regions[ridx]->m_growing_plus_tip_list.begin();
        while (tip_number > 0)
        {
            tip_tag++;
            tip_number--;
        }
    }
    else
    {
        tip_tag = --(m_geometry->m_regions[ridx]->m_growing_plus_tip_list.end());
        while (tip_number < m_geometry->m_regions[ridx]->m_growing_plus_tip_list.size() - 1)
        {
            tip_tag--;
            tip_number++;
        }
    }

    (*tip_tag)->m_mt->catastrophe();

    return;
}

void System::handle_rescue_event()
{
    // TODO
    // Fix the logic below and the types of these variables
    // to eliminate the need for static_cast.
    size_t ridx(0);
    unsigned long int tip_number(0);

    RegionMTTipTag tip_tag;

    // select the number of the tip to be rescued
    tip_number = m_random_gen.rand_int(m_shrinking_mts.size() - 1);

    // determine the region this tip is located : binary search algorithm
    if (tip_number < static_cast<unsigned long long>(m_shrinking_mts.size() / 2))
    {
        while (tip_number >= m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size())
        {
            tip_number -= m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size();
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tip_number = m_shrinking_mts.size() - 1 - tip_number;
        ridx = m_geometry->m_regions.size() - 1;
        while (tip_number >= m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size())
        {
            tip_number -= m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size();
            ridx--;
        }
        // flip the direction of the result for the next step
        tip_number = m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size() - 1 - tip_number;
    }

    // determine which tip in the selected region : binary search algorithm
    if (tip_number < m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size() / 2)
    {
        tip_tag = m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.begin();
        while (tip_number > 0)
        {
            tip_tag++;
            tip_number--;
        }
    }
    else
    {
        tip_tag = --(m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.end());
        while (tip_number < m_geometry->m_regions[ridx]->m_shrinking_plus_tip_list.size() - 1)
        {
            tip_tag--;
            tip_number++;
        }
    }

    // rescue the selected tip
    (*tip_tag)->m_mt->rescue();

    return;
}

void System::handle_severing_event()
{
    double cut_length(0.0);
    Segment* cut_seg(nullptr);
    random_position_on_microtubule(cut_length, cut_seg);

    // sever the MT now !
    if (cut_seg != nullptr)
    {
        cut_seg->m_mt->sever(cut_seg, cut_length);
    }
    else
    {
        cerr << "ERROR (non-fatal): severing position outside total MT length. Ignoring event.\n";
    }

    return;
}

void System::random_position_on_microtubule(double& random_pos, Segment*& random_seg)
{
    // TODO
    // Fix the logic below and the type of this variable
    // to eliminate the need for static_cast.
    long long ridx(0);
    // select position at which severing will take place
    double cut_length = m_random_gen.rand_dbl_exc(m_total_length);
    // first, select relevant region. Go from both sides to reduce time.
    if (cut_length < 0.5 * m_total_length)
    {
        ridx = 0;
        while (true)
        {
            m_geometry->m_regions[ridx]->update_region_length();
            if (cut_length < m_geometry->m_regions[ridx]->m_total_length)
            {
                break;
            }
            cut_length -= m_geometry->m_regions[ridx]->m_total_length;
            ridx++;
            if (ridx == static_cast<long long>(m_geometry->m_regions.size()))
            {
                // apparently, cutLength was *just* too long...
                cut_length -= ZERO_CUTOFF;
#ifdef DBG_ACID_TEST
                if (cutLength > geometry->regions[ridx - 1]->totalLength)
                {
                    cerr << "ERROR: cannot locate intersection location\n";
                    exit(-1);
                }
#endif
                break;
            }
        }
    }
    else
    {
        ridx = static_cast<long long>(m_geometry->m_regions.size() - 1);
        cut_length = m_total_length - cut_length;
        while (true)
        {
            m_geometry->m_regions[ridx]->update_region_length();
            if (cut_length < m_geometry->m_regions[ridx]->m_total_length)
            {
                break;
            }
            cut_length -= m_geometry->m_regions[ridx]->m_total_length;
            ridx--;
            if (ridx == -1)
            {
                // apparently, cutLength was *just* too long...
                cut_length -= ZERO_CUTOFF;

#ifdef DBG_ACID_TEST
                if (cutLength > geometry->regions[0]->totalLength)
                {
                    cerr << "ERROR: cannot locate intersection location\n";
                    exit(-1);
                }
#endif
                break;
            }
        }
        // and wrap it back
        cut_length = m_geometry->m_regions[ridx]->m_total_length - cut_length;
    }

    // Now that the area has been selected, cycle over trajectories and their segments
    Trajectory* trptr(nullptr);
    list<Segment*>::iterator seg;
    double temp(0.0);
    bool stop_flag(false);
    // Again, let the search direction depend on the expected position
    if (cut_length < 0.5 * (m_geometry->m_regions[ridx]->m_total_length))
    {
        // search in the forward direction. Start with the first trajectory
        trptr = m_geometry->m_regions[ridx]->m_trajectories.first();
        while (true)
        {
            // get the first segment on this trajectory
            seg = trptr->m_segments.begin();
            // and loop over all segments
            while (seg != trptr->m_segments.end())
            {
                // make sure the lengths are accurate
                (*seg)->m_mt->update_length();
                // break when the correct segment is found
                if (cut_length < (temp = (*seg)->length()))
                {
                    stop_flag = true;
                    break;
                }
                else
                {
                    // if not, continue
                    cut_length -= temp;
                    seg++;
                }
            }

            if (stop_flag)
            {
                break;
            }

            // if we ran out of trajectories, something is wrong. This could happen due to numerical inaccuracies.
            if (trptr->next() == nullptr)
            {
                cerr << "Rare event: cutting on the edge\n";
                seg--;
                cut_length -= ZERO_CUTOFF;
                break;
            }
            trptr = trptr->next();
        }
    }
    else
    {
        // flip the length and trajectory search orders
        cut_length = m_geometry->m_regions[ridx]->m_total_length - cut_length;
        trptr = m_geometry->m_regions[ridx]->m_trajectories.last();
        while (true)
        {
            // also, start with the last segment (one step beyond, in this case)
            seg = (trptr->m_segments.end());
            do
            {
                seg--;
                // make sure that the segment length is correct
                (*seg)->m_mt->update_length();
                // break if the requested position lies on this segment
                if (cut_length < (temp = (**seg).length()))
                {
                    stop_flag = true;
                    break;
                }
                else
                {
                    cut_length -= temp;
                }
            } while (seg != trptr->m_segments.begin());

            if (stop_flag)
            {
                break;
            }

            if (trptr->previous() == nullptr)
            {
                cerr << "Rare event: cutting on the edge\n";
                seg++;
                cut_length -= ZERO_CUTOFF;
                break;
            }
            trptr = trptr->previous();
        }
        cut_length = (*seg)->length() - cut_length;
    }

    random_seg = *seg;
    random_pos = cut_length;
    return;
}

#ifdef CROSS_SEV

// TODO: verify this function
void System::handle_severing_at_cross_event()
{
    Microtubule* cut_mt(nullptr);
    IntersectionItr cut_is;
    OccupiedIntersection* cut_occ_is(nullptr);
    Segment* cut_seg(nullptr);
    int cut_seg_number(0);
    list<Segment*>::iterator test_seg;
    list<Segment*>::iterator final_test_seg;

    cut_occ_is = m_occupied_intersection_list.element_address(m_random_gen.rand_int(m_occupied_intersection_list.size() - 1));
    if (m_p.m_cross_severing_top || (m_random_gen.rand_int(1) == 0))
    {
        // if MT (bundle) on top should be cut, select that one, or (if random) pick it with 50% chance.
        // Note that no effort is made to ensure proportional selection for bundles of different thickness.
        cut_is = cut_occ_is->m_intersection_to_cut;
    }
    else
    {
        cut_is = cut_occ_is->m_intersection_to_cut->second.m_mirror;
    }
    double diff_angle = cut_is->second.m_other_trajectory->m_base.m_region->intersection_angle(
    cut_is->second.m_other_trajectory, cut_is->second.m_mirror->second.m_other_trajectory);
    if (diff_angle > 0.5 * PI)
    {
        diff_angle = PI - diff_angle;
    }
    // don't do anything if the angle between MTs is smaller than the cutoff angle
    if (diff_angle <= (m_p.m_cross_severing_start_angle) * PI / 180.)
    {
        return;
    }

    /*
    IMPORTANT NOTE: DO NOT REMOVE OCCUPIED INTERSECTION HERE! IT WILL BE REMOVED AUTOMATICALLY IN A LATER EVENT
    cutIS->second.occupancy--;		// do not decrease: will happen very soon in a new event
    {
        removeOccupiedIntersection(cutIS->second);
    }
    */
    // find the segment to be cut
    cut_seg_number = m_random_gen.rand_int(cut_is->second.m_occupancy - 1); // nummers 0 - occupancy -1 .
    if (cut_seg_number < cut_is->second.m_occupancy / 2)                // search from the beginning.
                                                    // In principle, all searches could start from the beginning. Also
                                                    // starting from the end should save (on average) half of the time
    {
        test_seg = cut_is->second.m_mirror->second.m_other_trajectory->m_segments.begin();
        final_test_seg = cut_is->second.m_mirror->second.m_other_trajectory->m_segments.end();
        cut_seg_number++;
        while (test_seg != final_test_seg) // NOTE: test should not be necessary "while(true)"
        {
            if ((**test_seg).crosses_intersection(cut_is))
            {
                cut_seg_number--;
                if (cut_seg_number == 0)
                {
                    cut_seg = *test_seg;
                    break;
                }
            }
            test_seg++;
        }
    }
    else // search from the end
    {
        test_seg = cut_is->second.m_mirror->second.m_other_trajectory->m_segments.end();
        final_test_seg = cut_is->second.m_mirror->second.m_other_trajectory->m_segments.begin();
        while (test_seg != final_test_seg)
        {
            test_seg--;
            if ((**test_seg).crosses_intersection(cut_is))
            {
                cut_seg_number++;
                if (cut_seg_number == cut_is->second.m_occupancy)
                {
                    cut_seg = *test_seg;
                    break;
                }
            }
        }
    }
    if (cut_seg != nullptr)
    {
        cut_mt = cut_seg->m_mt;
        cut_mt->update_length();
        if (cut_mt != nullptr)
        {
            cut_mt->sever_at_cross(cut_is, cut_seg);
        }
        else
        {
            cerr << "ERROR (non-fatal): severing position outside total MT length. Ignoring event.\n";
        }
    }
    else
    {
        cerr << "ERROR (non-fatal): Did not find the right segment. Ignoring event.\n";
    }

    return;
}
#endif

void System::handle_nucleation_event()
{
    // if there is a finite tubulin pool and vPlus is essentially the same as vTM : leading to moving microtubules with
    // zero length, causing ordering problems on collisions
    if ((m_p.m_restricted_pool != 0) && (m_p.m_v_tm / (m_p.m_v_plus * (1. - m_total_length / (m_p.m_pool_density * m_geometry->m_area))) > 0.9999))
    {
        return;
    }

    // no nucleation at forbidden zone
    TrajectoryVector tv;

    // get random nucleation point (vector)
    SurfaceVector sv = m_geometry->random_surface_vector();

    // forbiden zone and parameter modified : no microtubule created
    if (sv.m_region->m_face_tag == -1)
    {
        return;
    }

    double override_angle(0.0);

    switch (m_p.m_nucleation_type)
    {
            // select isotropic nucleation
        case NUC_ISOTROPIC:
            override_angle = sv.m_angle;
            break;
        default:
            break;
    }

    sv.m_angle = override_angle;
    tv = m_geometry->create_trajectory(sv);
    m_growing_mts.create(this, tv);

    return;
}

void System::determine_stochastic_event()
{
    // for(regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);

    double next_event_interval(0.0);
    double inv_total_rate(0.0);
    double type_selector(0.0);
    double log_random(0.0);
    double temp(0.0);

    if ((m_growing_mts.size() == 0) || (m_p.m_restricted_pool == 0))
    {
        double const_rate(0.0);
        double slope_rate(0.0);

        const_rate = m_p.m_k_nuc * m_geometry->m_area + m_p.m_k_res * m_shrinking_mts.size() + m_p.m_k_sev * m_total_length;

        for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
        {
            const_rate += m_p.m_k_cat * m_p.m_region_kcat_multiplier[d_tag] * m_growing_tips_reg[d_tag];
        }

#ifdef CROSS_SEV
        const_rate += m_p.m_k_cross * m_occupied_intersection_list.size();
#endif

        slope_rate = m_p.m_k_sev * ((m_p.m_v_plus - m_p.m_v_tm) * m_growing_mts.size() + (m_p.m_v_min - m_p.m_v_tm) * m_shrinking_mts.size());

        // generate a random number larger than 0 with an exponential distribution
        log_random = -log(m_random_gen.rand_dbl_exc());

        if (abs(slope_rate) < ZERO_CUTOFF)
        {
            next_event_interval = log_random / const_rate;
        }
        else
        {
            if ((temp = pow(const_rate, 2) + 2 * slope_rate * log_random) < 0)
            {

                // choose a time that is always larger than the disappearance time of all MTs
                next_event_interval = VERY_LARGE;
            }

            else
            {
                next_event_interval = (-const_rate + sqrt(temp)) / slope_rate;
            }
        }

        // it will be negative for the case with no stochastic events
        inv_total_rate = 1.0 / (const_rate + slope_rate * next_event_interval);
    }

    else
    {
        double inv_l_max = 1.0 / (m_p.m_pool_density * m_geometry->m_area);
        double alpha = (m_p.m_v_plus - m_p.m_v_tm) * m_growing_mts.size() + (m_p.m_v_min - m_p.m_v_tm) * m_shrinking_mts.size();
        double beta = m_p.m_v_plus * m_growing_mts.size() * inv_l_max;
        double l_inf = alpha / beta;

        double base_rate = m_p.m_k_nuc * m_geometry->m_area + m_p.m_k_res * m_shrinking_mts.size();
        for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
        {
            base_rate += m_p.m_k_cat * m_p.m_region_kcat_multiplier[d_tag] * m_growing_tips_reg[d_tag];
        }

#ifdef CROSS_SEV
        base_rate += m_p.m_k_cross * m_occupied_intersection_list.size();
#endif

        double x = m_p.m_k_sev * (m_total_length - l_inf) / (base_rate + m_p.m_k_sev * l_inf);
        double y = -beta * log(m_random_gen.rand_dbl_exc()) / (base_rate + m_p.m_k_sev * l_inf);

        next_event_interval = (-x + y + lambert_w(x * exp(x - y))) / beta;

        double total_rate = base_rate + m_p.m_k_sev * (l_inf + (m_total_length - l_inf) * exp(-beta * next_event_interval));
        inv_total_rate = 1.0 / total_rate;
    }

#ifdef DBG_ACID_TEST
    if ((nextEventInterval < 0) || (nextEventInterval != nextEventInterval))
    {
        cerr << "Invalid deterministic event time calculated. Exiting\n";
        exit(-1);
    }
#endif

    m_next_stochastic_event_time = m_system_time + next_event_interval;

    // even if  invTotalRate is negative! due to the timing, the event will never be executed.
    type_selector = m_random_gen();

    double type_select_reg_source(0.0);
    for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
    {
        type_select_reg_source += m_p.m_k_cat * m_p.m_region_kcat_multiplier[d_tag] * m_growing_tips_reg[d_tag];
    }

    if (type_selector < m_p.m_k_nuc * m_geometry->m_area * inv_total_rate)
    {
        m_next_stochastic_event_type = NUCLEATION;
    }

    // select catastrophe event
    else if (type_selector < (m_p.m_k_nuc * m_geometry->m_area + type_select_reg_source) * inv_total_rate)
    {
        m_next_stochastic_event_type = CATASTROPHE;
    }

    // select rescue event
    else if (type_selector
             < (m_p.m_k_nuc * m_geometry->m_area + type_select_reg_source + m_p.m_k_res * m_shrinking_mts.size()) * inv_total_rate)
    {
        m_next_stochastic_event_type = RESCUE;
    }

// select severing event
#ifdef CROSS_SEV
    else if (type_selector < (m_p.m_k_nuc * m_geometry->m_area + type_select_reg_source + m_p.m_k_res * m_shrinking_mts.size()
                             + m_p.m_k_cross * m_occupied_intersection_list.size())
                            * inv_total_rate)
    {
        m_next_stochastic_event_type = SEVERING_AT_CROSS;
    }
#endif

    // select katanin induced event
    else
    {
        m_next_stochastic_event_type = KATANIN;
    }

#ifdef DBG_EVENT
    cerr << "determined next stochastic event. type=" << nextStochasticEventType << ", time=" << nextStochasticEventTime
         << "\n";
#endif

    return;
}

void System::handle_global_event(DeterministicEvent& event)
{
    // better to be updated
    update_all();

    switch (event.m_global_type)
    {

        case MEASURE:
            break;
        case STATUS:
            if (!integrity_check())
            {
                cerr << "Integrity check failed. Exiting.\n";
                exit(-1);
            }

            perform_measurement();
            m_next_status_event_time = m_system_time + m_p.m_measurement_interval;
            m_time_queue.push_global(m_next_status_event_time, STATUS);
            break;

        case SNAPSHOT:
            // declaring variable within case: so the local blocking
            {
                OrderParameters order;

                m_geometry->get_order_parameters(order);

                // write a movie frame
                m_movie_file << order.m_rdirector[0] << "," << order.m_rdirector[1] << "," << order.m_rdirector[2] << "#";
                m_movie_file << m_geometry->m_object_cm(0, 0) << "," << m_geometry->m_object_cm(1, 0) << ","
                          << m_geometry->m_object_cm(2, 0) << "#";
                m_movie_file << m_geometry->m_object_pa(0, 0) << "," << m_geometry->m_object_pa(1, 0) << ","
                          << m_geometry->m_object_pa(2, 0) << "*";

                m_geometry->output_snapshot(m_movie_file);

                // write a local order parameter heatmap frame
                m_geometry->output_order_heat_map(m_heat_map_file, order.m_local_order, order.m_sv);

                m_next_snapshot_event_time = m_system_time + m_p.m_movie_frame_interval;
                m_time_queue.push_global(m_next_snapshot_event_time, SNAPSHOT);
            }
            break;

        case STOP:
        {
            OrderParameters order;
            m_geometry->get_order_parameters(order);

            // division plane: polygon edge-list & area
            Eigen::Vector3d normal_ppb(order.m_rdirector[0], order.m_rdirector[1], order.m_rdirector[2]);
            m_geometry->m_area_ppb = intersecting_polygon(m_geometry->m_ppb,
                                                    m_geometry->m_global_vertex,
                                                    m_geometry->m_regions,
                                                    m_geometry->m_ppb_edge_list,
                                                    normal_ppb,
                                                    m_geometry->m_nucleous_position,
                                                    m_p.m_geometry);

            if (m_p.m_show_output == 1)
            {
                cout << "Division plane area: " << m_geometry->m_area_ppb << endl;
            }

            m_stop_signal = true;
        }
        break;

        case PARAMETER_CHANGE:
        {
            cout << "Changing parameters. Reading from file [" << m_p.m_new_parameter_file << "]\n";
            flush_and_reload(false);
            string old_dir = m_p.m_output_dir;
            if (m_p.reinitialize(m_p.m_new_parameter_file.c_str()))
            {
                if (m_p.m_output_dir != old_dir)
                {
                    cout << "Switching directories.\n";
                    write_measurements_to_file(0);
                    close_files();
                    initialize_output();
                }
                m_p.write_to_file();

                if (m_p.m_movie_enabled)
                {
                    if (m_next_snapshot_event_time < ZERO_CUTOFF)
                    {
                        m_time_queue.push_global(m_p.m_movie_frame_interval, SNAPSHOT);
                        m_next_snapshot_event_time = m_p.m_movie_frame_interval;
                    }
                }
                else
                {
                    m_next_snapshot_event_time = 0;
                }

                if (m_p.m_new_parameter_read_interval > ZERO_CUTOFF)
                {
                    m_time_queue.push_global(m_p.m_new_parameter_read_interval, PARAMETER_CHANGE);
                    m_next_parameter_event_time = m_p.m_new_parameter_read_interval;
                }
            }
            else
            {
                cerr << "Failed to read parameters from file. Exiting.\n";
                emergency_break();
                exit(-1);
            }

            m_time_queue.push_global(m_p.m_new_parameter_read_interval, PARAMETER_CHANGE);
            m_next_parameter_event_time = m_p.m_new_parameter_read_interval;
            break;
        }

        case SIGNAL_PPB:
            cout << " Starting formation of PPB ...\n";
            m_p.write_to_file();
            m_p.m_pp_bk_nuc_fraction = 0.8;
            cout << m_p.m_pp_bk_nuc_fraction << " : " << m_p.m_ppb << endl;
            OrderParameters order;
            m_geometry->get_order_parameters(order);
            Eigen::Vector3d normal_ppb(order.m_rdirector[0], order.m_rdirector[1], order.m_rdirector[2]);
            m_geometry->m_area_ppb = intersecting_polygon(m_geometry->m_ppb,
                                                    m_geometry->m_global_vertex,
                                                    m_geometry->m_regions,
                                                    m_geometry->m_ppb_edge_list,
                                                    normal_ppb,
                                                    m_geometry->m_nucleous_position,
                                                    m_p.m_geometry);

            for (int d_tag = 0; d_tag <= m_p.m_face_number; d_tag++)
            {
                m_growing_tips_reg[d_tag] = 0;
            }

            for (size_t eno = 0; eno < m_geometry->m_regions.size(); eno++)
            {
                if (m_geometry->m_regions[eno]->m_poly_intersect_mark == 1)
                {
                    m_geometry->m_regions[eno]->m_face_tag = 1;

                    for (size_t tip = 0; tip < m_geometry->m_regions[eno]->m_growing_plus_tip_list.size(); tip++)
                    {
                        m_growing_tips_reg[1] += 1;
                    }

                    m_geometry->m_patch_area[1] += m_geometry->m_regions[eno]->m_area;
                    m_geometry->m_regions_index[1].push_back(eno);
                }

                else
                {
                    m_geometry->m_regions[eno]->m_face_tag = 2;

                    for (size_t tip = 0; tip < m_geometry->m_regions[eno]->m_growing_plus_tip_list.size(); tip++)
                    {
                        m_growing_tips_reg[2] += 1;
                    }

                    m_geometry->m_patch_area[2] += m_geometry->m_regions[eno]->m_area;
                    m_geometry->m_regions_index[2].push_back(eno);
                }
            }
            m_time_queue.push_global(VERY_LARGE, SIGNAL_PPB);
            m_next_ppb_event_time = VERY_LARGE;

            break;
    }

    return;
}

EventDescriptorIndex System::register_event_descriptor(EventDescriptor* ei)
{
    // maybe use a pool of numbers instead of an increasing number....this could wrap around
    EventDescriptorIndex index = m_event_descriptor_id++;
    m_event_descriptor_map[index] = ei;

    if (m_event_descriptor_id == numeric_limits<EventDescriptorIndex>::max() - 10)
    {
        m_time_queue.push_global(m_system_time, STOP);
    }

    return index;
}

void System::unregister_event_descriptor(EventDescriptorIndex ei)
{
    m_event_descriptor_map.erase(m_event_descriptor_map.find(ei));
    return;
}

EventDescriptor* System::get_event_descriptor(EventDescriptorIndex idx)
{
    map<EventDescriptorIndex, EventDescriptor*>::iterator itr = m_event_descriptor_map.find(idx);
    if (itr == m_event_descriptor_map.end())
    {
        return nullptr;
    }
    else
    {
        return itr->second;
    }
}

#ifdef CROSS_SEV
void System::remove_occupied_intersection(Intersection& is)
{
    m_occupied_intersection_list.remove_element(is.m_occupied_list_ptr->m_index);
    is.m_occupied_list_ptr = nullptr;
    is.m_mirror->second.m_occupied_list_ptr = nullptr;
    return;
}

void System::add_occupied_intersection(IntersectionItr is)
{
    is->second.m_occupied_list_ptr = m_occupied_intersection_list.create(is);
    is->second.m_mirror->second.m_occupied_list_ptr = is->second.m_occupied_list_ptr;
    return;
}
#endif

void System::make_binomial_table(void)
{
    int i(0), j(0);
    double temp(0.0);
    double factorial[MAXBINOM];
    temp = 1;
    factorial[0] = 1;
    for (i = 1; i < MAXBINOM; i++)
    {
        temp *= i;
        factorial[i] = temp;
    }
    for (i = 0; i < MAXBINOM; i++)
    {
        for (j = 0; j <= i; j++)
        {
            m_binomial_table[i][j] = factorial[i] / (factorial[j] * factorial[i - j]);
        }
    }
    return;
}
