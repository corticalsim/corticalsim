#include "mesh.hpp"
#include "region.hpp"
#include "system.hpp"
#include "mt_tip.hpp"
#include "parameters.hpp"

void System::collision_probabilities(double angle, double& pcat, double& pzip)
{
    switch (p.interaction_type)
    {
        case InteractionType::ZIP_FIRST:

            // angle < magic angle
            if (angle < PI * p.magic_angle / 180.0)
            {
                pcat = 0.0;

                // zippering enabled
                if (p.zippering_enabled)
                {
                    pzip = p.zip_fraction;
                }
                else
                {
                    pzip = 0.0;
                }
                break;
            }

            pzip = 0.0;

            // angle > p.catStartAngle and induced catastrophes enabled
            if (p.catastrophes_enabled && (angle > p.cat_start_angle / 180.0 * PI))
            {
                // proportional catastrophe enabled
                if (p.proportional_catastrophes)
                {
                    pcat = (angle - p.cat_start_angle / 180.0 * PI) / (PI / 2. - p.cat_start_angle / 180.0 * PI)
                           * p.induced_catastrophe_fraction;
                }

                // proportional catastrophe disabled
                else
                {
                    pcat = p.induced_catastrophe_fraction;
                }
            }
            else
            {
                pcat = 0.0;
            }

            break;

        case InteractionType::CAT_FIRST:

            // angle > p.catStartAngle and induced catastrophes enabled
            if (p.catastrophes_enabled && (angle > p.cat_start_angle / 180.0 * PI))
            {
                // proportional catastrophe enabled
                if (p.proportional_catastrophes)
                {
                    pcat = (angle - p.cat_start_angle / 180.0 * PI) / (PI / 2. - p.cat_start_angle / 180.0 * PI)
                           * p.induced_catastrophe_fraction;
                }

                // proportional catastrophe disabled
                else
                {
                    pcat = p.induced_catastrophe_fraction;
                }
            }
            else
            {
                pcat = 0.0;
            }

            // angle < magic angle and zippering enabled
            if (p.zippering_enabled && (angle < p.magic_angle / 180.0 * PI))
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
            total_sum += (npar + 1 - i) * binomial_table[ncoll][i - j] * binomial_table[max(0, i - 1)][j]
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
                total_sum += binomial_table[j][i - k] * binomial_table[i][k] * pow(x_single, 2 * k + j - i)
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
    double rand = random_gen.rand_dbl_exc();

    // get the available Pcat and Pzip
    collision_probabilities(angle, pcat, pzip);

    switch (p.bundle_type)
    {
        // TODO
        // This needs to be fixed properly.
        // Adding this case or `default` at the end causes an error becasue
        // of the initialisation of multiPx & multiPz in the last case.
        case BundleType::COUNT_LAST:
            break;
        case BundleType::NO_ZIP:

            if (n_parallel > 0)
            {
                if (rand < pcat)
                {
                    return CollisionType::INDUCED_CATASTROPHE;
                }
            }
            // break omitted on purpose!

        case BundleType::STICKY:

            if (n_parallel > 0)
            {
                return CollisionType::CROSSOVER;
            }
            // break omitted on purpose!

        case BundleType::SIMPLE:

            // induced catastrophe selected
            if (rand < pcat)
            {
                return CollisionType::INDUCED_CATASTROPHE;
            }

            // zippering selected
            else if (rand < pcat + pzip)
            {
                return CollisionType::ZIPPER;
            }

            // cross over slected
            else
            {
                return CollisionType::CROSSOVER;
            }

            // though never reached
            // break;

        case BundleType::NCOLLISION:

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
                return CollisionType::INDUCED_CATASTROPHE;
            }

            else if (rand < pcat + pzip)
            {
                return CollisionType::ZIPPER;
            }

            else
            {
                return CollisionType::CROSSOVER;
            }

            // though never reached
            // break;

        case BundleType::MULTI_COLLISION:

            double pcross = 1.0 - pcat - pzip;
            double multi_px = multi_pcross(min(n_parallel, MAXBINOM - 1), min(n_cross, MAXBINOM - 1), pcross, pzip);
            double multi_pz = multi_pzip(min(n_parallel, MAXBINOM - 1), min(n_cross, MAXBINOM - 1), pcross, pzip);

            if (rand < multi_px)
            {
                return CollisionType::CROSSOVER;
            }

            else if (rand < multi_px + multi_pz)
            {
                return CollisionType::ZIPPER;
            }

            else
            {
                return CollisionType::INDUCED_CATASTROPHE;
            }

            // though never reached
            // break;
    }

    // should never occur
    return CollisionType::CROSSOVER;
}

System::System(char* par_file):
    geometry(nullptr),
    p(this),
    system_time(0),
    system_time_offset(0),
    current_time_tag(0),
    total_length(0),
    mem_usage(0),
    stop_signal(false),
    wall_clock_start_time(time(nullptr)),
    total_s_event_count(0),
    total_valid_d_event_count(0),
    total_invalid_d_event_count(0),
    total_length_severing_count(0),
    total_intersection_severing_count(0),
    boundary_crossing_count(0),
    total_zipper_count(0),
    total_crossover_count(0),
    total_induced_catastrophe_count(0),
    count_segments(0),
    count_trajectories(0),
    count_intersections(0),
    time_queue(this, &System::identity, &System::identity),
    v_plus_queue(this, &System::v_plus_to_time, &System::time_to_v_plus),
    next_status_event_time(0),
    next_snapshot_event_time(0),
    event_descriptor_id(0),
    event_id(0),
    data_saved(false)
{
    /// read parameters
    p.initialize(par_file);

    // initialize number of growing tips in different faces
    for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
    {
        growing_tips_reg.push_back(0);
    }

    // load geometry
    geometry = new TriMeshGeometry(this);

    /// initialize the random seed with the desired one
    random_gen.seed(p.seed);

    /// used in multiPcross & multiPzip events
    make_binomial_table();

    /// perform further initialization
    initialize_output();
    p.write_to_file();

    /// start creating MT events
    determine_stochastic_event();
    time_queue.push_global(p.measurement_interval, GlobalEventType::STATUS);
    next_status_event_time = p.measurement_interval;

    /// if snapShot events are requested
    if (p.movie_enabled)
    {
        time_queue.push_global(p.movie_frame_interval, GlobalEventType::SNAPSHOT);
        next_snapshot_event_time = p.movie_frame_interval;
    }

    /// if parameter file is updated on run
    if (p.ppb)
    {
        time_queue.push_global(p.pp_bforming_time, GlobalEventType::SIGNAL_PPB);
        next_ppb_event_time = p.pp_bforming_time;
    }

    /// if parameter file is updated on run
    if (p.new_parameter_read_interval > ZERO_CUTOFF)
    {
        time_queue.push_global(p.new_parameter_read_interval, GlobalEventType::PARAMETER_CHANGE);
        next_parameter_event_time = p.new_parameter_read_interval;
    }

    return;
}

System::~System()
{
    // destroy all MTs
    growing_mts.remove_all();
    shrinking_mts.remove_all();

    return;
}

void System::run(double sim_time, string& c_s, string& m_s)
{
    int wall_time_polling_counter(0);
    int memory_polling_counter(0);
    time_queue.push_global(sim_time - system_time_offset, GlobalEventType::STOP);

    if (p.show_output == 1)
    {
        cout << "===============================================" << endl;
        cout << " Time stamped simulation output " << endl;
        cout << "===============================================" << endl;
    }

    while (!stop_signal)
    {
        if (++wall_time_polling_counter == CLOCK_POLLING_INTERVAL)
        {
            if (difftime(time(nullptr), wall_clock_start_time) >= p.wall_clock_limit)
            {
                // time's up; send a stop signal
                cerr << "Wall clock limit (" << p.wall_clock_limit << " seconds) reached. Stopping.\n";
                time_queue.push_global(system_time, GlobalEventType::STOP);

                // write final snapshot when running out of time.
                if (p.movie_enabled)
                {
                    geometry->output_snapshot(movie_file);
                }

                // emergency exit the program
                emergency_break();
                exit(-1);
            }
            wall_time_polling_counter = 0;
        }

        if (++memory_polling_counter == MEMORY_POLLING_INTERVAL)
        {
            if (estimate_memory_footprint() > 1024 * 1024 * p.memory_limit)
            {
                cerr << "Memory limit (approx. " << p.memory_limit << " MB) reached. Stopping.\n";
                time_queue.push_global(system_time, GlobalEventType::STOP);

                // emergency exit the program
                emergency_break();
                exit(-1);
            }
            memory_polling_counter = 0;
        }

        if (event_id > QUEUE_FLUSH_INTERVAL)
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
    total += (growing_mts.size() + shrinking_mts.size()) * (sizeof(Microtubule) + 2 * sizeof(MTTip*));
    total += event_descriptor_map.size() * (sizeof(pair<EventDescriptorIndex, EventDescriptor*>) + 3 * sizeof(void*));
    total += count_segments * (sizeof(Segment) + sizeof(Segment*));
    total += count_trajectories * (sizeof(Trajectory));
    total += count_intersections * (2 * (sizeof(pair<double, Intersection>) + 3 * sizeof(void*)));
    total += (time_queue.queue.size() + v_plus_queue.queue.size()) * (sizeof(DeterministicEvent) + 2 * sizeof(void*));
    total += measurement_history.size() * sizeof(Measurement);

    return total;
}

bool System::integrity_check()
{
    bool valid(true);

    // integrity check of the gemoetry embeded components
    if (!geometry->integrity_check())
    {
        valid = false;
    }

    // integrity check of each growing MTs
    Microtubule* mt = growing_mts.first();
    while (mt != nullptr)
    {
        if (!mt->integrity_check())
        {
            valid = false;
        }

        mt = mt->next();
    }

    // integrity check of each shrinking MTs
    mt = shrinking_mts.first();
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
    double time_offset(system_time);

    // update all microtubule lengths
    update_all();

    // recalculate total lengths to avoid drift (in DBG mode, also check for large deviations)
    double region_length(0.);
    double system_length(0);
    Trajectory* tptr(nullptr);

    for (size_t ridx = 0; ridx < geometry->regions.size(); ridx++)
    {
        region_length = 0;
        tptr = geometry->regions[ridx]->trajectories.first();

        while (tptr != nullptr)
        {
            region_length += tptr->segment_length();
            tptr = tptr->next();
        }

        if (abs(geometry->regions[ridx]->total_length - region_length) > ZERO_CUTOFF * max(100., region_length))
        {
            cerr << "Unacceptable drift in region length (drift = "
                 << abs(geometry->regions[ridx]->total_length - region_length)
                 << "). Exiting. [try lowering flush interval]\n";
            emergency_break();
            exit(-1);
        }
        geometry->regions[ridx]->total_length = region_length;
        system_length += region_length;
    }

    if (abs(total_length - system_length) > ZERO_CUTOFF * max(1., total_length))
    {
        cerr << "Unacceptable drift in system length (drift = " << abs(total_length - system_length)
             << "). Exiting. [try lowering flush interval]\n";
        emergency_break();
        exit(-1);
    }

    total_length = system_length;

    // reset master timer
    system_time_offset += time_offset;

    system_time = 0;
    time_queue.flush();
    v_plus_queue.flush();
    event_descriptor_map.clear();
    event_id = 0;
    event_descriptor_id = 0;

    // now, reload all events (this is a PROBLEM for global events that have not been saved. They need to be re-inserted
    // manually)
    determine_stochastic_event();
    time_queue.push_global(p.stop_time - system_time_offset, GlobalEventType::STOP);
    next_status_event_time -= time_offset;
    time_queue.push_global(next_status_event_time, GlobalEventType::STATUS);

    if (refresh_parameter_event)
    {
        next_parameter_event_time -= time_offset;
        time_queue.push_global(next_parameter_event_time, GlobalEventType::PARAMETER_CHANGE);
    }

    if (p.movie_enabled)
    {
        next_snapshot_event_time -= time_offset;
        time_queue.push_global(next_snapshot_event_time, GlobalEventType::SNAPSHOT);
    }

    if (p.ppb)
    {
        next_ppb_event_time -= time_offset;
        time_queue.push_global(next_ppb_event_time, GlobalEventType::SIGNAL_PPB);
    }

    Microtubule* mt(nullptr);
    Segment* seg(nullptr);
    mt = growing_mts.first();

    while (mt != nullptr)
    {
        // this can go into a change time function
        mt->nucleation_time -= time_offset;
        seg = mt->segments.first();

        while (seg != nullptr)
        {
            seg->nucleation_time -= time_offset;
            seg = seg->next();
        }

        mt->plus.event.index = register_event_descriptor(&(mt->plus.event));
        mt->minus.event.index = register_event_descriptor(&(mt->minus.event));
        mt->disappear_event.index = register_event_descriptor(&(mt->disappear_event));
        mt->plus.determine_event();
        mt->minus.determine_event();
        mt->set_disappear_event();
        mt = mt->next();
    }

    mt = shrinking_mts.first();

    while (mt != nullptr)
    {
        mt->nucleation_time -= time_offset;
        seg = mt->segments.first();

        {
            seg->nucleation_time -= time_offset;
            seg = seg->next();
        }

        mt->plus.event.index = register_event_descriptor(&(mt->plus.event));
        mt->minus.event.index = register_event_descriptor(&(mt->minus.event));
        mt->disappear_event.index = register_event_descriptor(&(mt->disappear_event));
        mt->plus.determine_event();
        mt->minus.determine_event();
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
    if (p.restricted_pool == 0)
    {
        return dist;
    }

    // for finite tubulin pool
    double lmax = (p.pool_density * geometry->area);
    double inv_l_max = 1.0 / (p.pool_density * geometry->area);

    double alpha = (p.v_plus - p.v_tm) * growing_mts.size() + (p.v_min - p.v_tm) * shrinking_mts.size();

    if (growing_mts.size() == 0)
    {
        if (shrinking_mts.size() == 0)
        {
            return dist / (1.0 - total_length * inv_l_max);
        }
        else
        {
            return ((lmax - total_length) - sqrt(pow(lmax - total_length, 2) - 2. * alpha * lmax * dist)) / alpha;
        }
    }

    double beta = p.v_plus * growing_mts.size() * inv_l_max;

    if ((shrinking_mts.size() == 0) && ((p.v_tm / p.v_plus) < ZERO_CUTOFF))
    {
        return -log(1.0 - growing_mts.size() * p.v_plus * dist / (lmax - total_length)) / beta;
    }

    double pv = (alpha / beta - total_length) / (lmax - alpha / beta);
    double qv = growing_mts.size() * p.v_plus * dist / (lmax - alpha / beta);

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
    if (p.restricted_pool == 0)
    {
        return time;
    }

    // for finite tubulin pool
    double inv_l_max = 1.0 / (p.pool_density * geometry->area);
    double alpha = (p.v_plus - p.v_tm) * growing_mts.size() + (p.v_min - p.v_tm) * shrinking_mts.size();

    if (growing_mts.size() == 0)
    {
        return ((1.0 - total_length * inv_l_max) - 0.5 * alpha * inv_l_max * time) * time;
    }

    double beta = p.v_plus * growing_mts.size() * inv_l_max;
    double lambda_inf = inv_l_max * alpha / beta;
    double lambda_0 = total_length * inv_l_max;

    return (1.0 - lambda_inf) * time - (lambda_0 - lambda_inf) * (1.0 - exp(-beta * time)) / beta;
}

void System::advance_time(double new_time)
{
#ifdef DBG_SYSTEM
    cout << "DBG/SYSTEM: advancing time by " << newTime - systemTime << " seconds. New system time is " << newTime
         << ".\n";
#endif

    // make a jump in time step (taking advantage of using event-driven algorithm)
    time_queue.advance_time(new_time);
    v_plus_queue.advance_time(new_time);
    system_time = new_time;
    total_length += time_queue.progression(current_time_tag)
                   * (p.v_min * shrinking_mts.size() - p.v_tm * (growing_mts.size() + shrinking_mts.size()))
                   + v_plus_queue.progression(current_time_tag) * (p.v_plus * growing_mts.size());

    if (++current_time_tag == POSITION_CACHE_SIZE)
    {
        current_time_tag = 0;
        update_all(true);
    }

    time_queue.store_time(current_time_tag);
    v_plus_queue.store_time(current_time_tag);

    return;
}

void System::update_all(bool force_update)
{
    // always nice to be updated in lengths of all regions and microtubules (convenient for snapshots, etc.)
    for (size_t ridx = 0; ridx < geometry->regions.size(); ridx++)
    {
        geometry->regions[ridx]->update_region_length(force_update);
    }

    Microtubule* mt;
    mt = growing_mts.first();

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

    mt = shrinking_mts.first();
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
    if (!time_queue.empty())
    {
        tq_value = time_queue.first_event_time();
    }

    if (!v_plus_queue.empty())
    {
        dq_value = v_plus_queue.first_event_time();
    }

    if ((next_stochastic_event_time < tq_value) && (next_stochastic_event_time < dq_value))
    {
        // execute stochastic event
        total_s_event_count++;
#ifdef DBG_EVENT
        cout << "DBG/Event: executing stochastic event [time=" << nextStochasticEventTime
             << ", type=" << nextStochasticEventType << "]\n";
#endif

        advance_time(next_stochastic_event_time);

        // execute next stochastic event
        switch (next_stochastic_event_type)
        {
            case StochasticEventType::CATASTROPHE:
                handle_catastrophe_event();
                break;
            case StochasticEventType::RESCUE:
                handle_rescue_event();
                break;
            case StochasticEventType::KATANIN:
                handle_severing_event();
                break;
#ifdef CROSS_SEV
            case StochasticEventType::SEVERING_AT_CROSS:
                handle_severing_at_cross_event();
                break;
#endif
            case StochasticEventType::NUCLEATION:
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
            event = time_queue.pop();
            next_event_time = tq_value;
        }
        else
        {
            event = v_plus_queue.pop();
            next_event_time = dq_value;
        }

        if (event.info_idx != -1)
        {
            map<EventDescriptorIndex, EventDescriptor*>::iterator event_itr = event_descriptor_map.find(event.info_idx);

            // whether the tip still exists and the tag still matches (valid event)
            if ((event_itr != event_descriptor_map.end()) && (event_itr->second->tag == event.tag))
            {
#ifdef DBG_EVENT
                cout << "DBG/Event: executing deterministic event [time=" << event.eventTimeDist
                     << ", type=" << eventItr->second->type << ", tag=" << event.tag << "]\n";
#endif

                total_valid_d_event_count++;
                advance_time(next_event_time);
                event_itr->second->mt->handle_event(event_itr->second);
                determine_stochastic_event();
            }

            // invalid event
            else
            {
                total_invalid_d_event_count++;
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

            advance_time(event.event_time_dist);
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
    for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
    {
        norm += growing_tips_reg[d_tag] * p.region_kcat_multiplier[d_tag];
    }

    double rand_gen = random_gen.rand();

    double norm_fraction(0.0);
    for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
    {
        norm_fraction += growing_tips_reg[d_tag] * p.region_kcat_multiplier[d_tag];
        if (rand_gen < (double)(norm_fraction) / norm)
        {
            domain_type = d_tag;
            total_relevant_tips = growing_tips_reg[d_tag];
            break;
        }
    }

    tip_number = random_gen.rand_int(total_relevant_tips - 1);

    // determine the region this tip is located in
    // depending on whether it's a high or low number, start scanning from the beginning or end
    if (tip_number < total_relevant_tips / 2)
    {
        while (true)
        {
            if (geometry->regions[ridx]->face_tag == domain_type)
            {
                if (tip_number < geometry->regions[ridx]->growing_plus_tip_list.size())
                {
                    break;
                }

                tip_number -= geometry->regions[ridx]->growing_plus_tip_list.size();
            }
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tip_number = total_relevant_tips - 1 - tip_number;
        ridx = geometry->regions.size() - 1;
        while (true)
        {
            if (geometry->regions[ridx]->face_tag == domain_type)
            {
                if (tip_number < geometry->regions[ridx]->growing_plus_tip_list.size())
                {
                    break;
                }
                tip_number -= geometry->regions[ridx]->growing_plus_tip_list.size();
            }
            ridx--;
        }
        // flip the direction of the result for the next step
        tip_number = geometry->regions[ridx]->growing_plus_tip_list.size() - 1 - tip_number;
    }

    // Now, select the tip with # 'tipNumber' from the current region.
    // Again, if the number is high within the number of tips in the region, start from the end
    if (tip_number < geometry->regions[ridx]->growing_plus_tip_list.size() / 2)
    {
        tip_tag = geometry->regions[ridx]->growing_plus_tip_list.begin();
        while (tip_number > 0)
        {
            tip_tag++;
            tip_number--;
        }
    }
    else
    {
        tip_tag = --(geometry->regions[ridx]->growing_plus_tip_list.end());
        while (tip_number < geometry->regions[ridx]->growing_plus_tip_list.size() - 1)
        {
            tip_tag--;
            tip_number++;
        }
    }

    (*tip_tag)->mt->catastrophe();

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
    tip_number = random_gen.rand_int(shrinking_mts.size() - 1);

    // determine the region this tip is located : binary search algorithm
    if (tip_number < static_cast<unsigned long long>(shrinking_mts.size() / 2))
    {
        while (tip_number >= geometry->regions[ridx]->shrinking_plus_tip_list.size())
        {
            tip_number -= geometry->regions[ridx]->shrinking_plus_tip_list.size();
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tip_number = shrinking_mts.size() - 1 - tip_number;
        ridx = geometry->regions.size() - 1;
        while (tip_number >= geometry->regions[ridx]->shrinking_plus_tip_list.size())
        {
            tip_number -= geometry->regions[ridx]->shrinking_plus_tip_list.size();
            ridx--;
        }
        // flip the direction of the result for the next step
        tip_number = geometry->regions[ridx]->shrinking_plus_tip_list.size() - 1 - tip_number;
    }

    // determine which tip in the selected region : binary search algorithm
    if (tip_number < geometry->regions[ridx]->shrinking_plus_tip_list.size() / 2)
    {
        tip_tag = geometry->regions[ridx]->shrinking_plus_tip_list.begin();
        while (tip_number > 0)
        {
            tip_tag++;
            tip_number--;
        }
    }
    else
    {
        tip_tag = --(geometry->regions[ridx]->shrinking_plus_tip_list.end());
        while (tip_number < geometry->regions[ridx]->shrinking_plus_tip_list.size() - 1)
        {
            tip_tag--;
            tip_number++;
        }
    }

    // rescue the selected tip
    (*tip_tag)->mt->rescue();

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
        cut_seg->mt->sever(cut_seg, cut_length);
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
    double cut_length = random_gen.rand_dbl_exc(total_length);
    // first, select relevant region. Go from both sides to reduce time.
    if (cut_length < 0.5 * total_length)
    {
        ridx = 0;
        while (true)
        {
            geometry->regions[ridx]->update_region_length();
            if (cut_length < geometry->regions[ridx]->total_length)
            {
                break;
            }
            cut_length -= geometry->regions[ridx]->total_length;
            ridx++;
            if (ridx == static_cast<long long>(geometry->regions.size()))
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
        ridx = static_cast<long long>(geometry->regions.size() - 1);
        cut_length = total_length - cut_length;
        while (true)
        {
            geometry->regions[ridx]->update_region_length();
            if (cut_length < geometry->regions[ridx]->total_length)
            {
                break;
            }
            cut_length -= geometry->regions[ridx]->total_length;
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
        cut_length = geometry->regions[ridx]->total_length - cut_length;
    }

    // Now that the area has been selected, cycle over trajectories and their segments
    Trajectory* trptr(nullptr);
    list<Segment*>::iterator seg;
    double temp(0.0);
    bool stop_flag(false);
    // Again, let the search direction depend on the expected position
    if (cut_length < 0.5 * (geometry->regions[ridx]->total_length))
    {
        // search in the forward direction. Start with the first trajectory
        trptr = geometry->regions[ridx]->trajectories.first();
        while (true)
        {
            // get the first segment on this trajectory
            seg = trptr->segments.begin();
            // and loop over all segments
            while (seg != trptr->segments.end())
            {
                // make sure the lengths are accurate
                (*seg)->mt->update_length();
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
        cut_length = geometry->regions[ridx]->total_length - cut_length;
        trptr = geometry->regions[ridx]->trajectories.last();
        while (true)
        {
            // also, start with the last segment (one step beyond, in this case)
            seg = (trptr->segments.end());
            do
            {
                seg--;
                // make sure that the segment length is correct
                (*seg)->mt->update_length();
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
            } while (seg != trptr->segments.begin());

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

    cut_occ_is = occupied_intersection_list.element_address(random_gen.rand_int(occupied_intersection_list.size() - 1));
    if (p.cross_severing_top || (random_gen.rand_int(1) == 0))
    {
        // if MT (bundle) on top should be cut, select that one, or (if random) pick it with 50% chance.
        // Note that no effort is made to ensure proportional selection for bundles of different thickness.
        cut_is = cut_occ_is->intersection_to_cut;
    }
    else
    {
        cut_is = cut_occ_is->intersection_to_cut->second.mirror;
    }
    double diff_angle = cut_is->second.other_trajectory->base.region->intersection_angle(
    cut_is->second.other_trajectory, cut_is->second.mirror->second.other_trajectory);
    if (diff_angle > 0.5 * PI)
    {
        diff_angle = PI - diff_angle;
    }
    // don't do anything if the angle between MTs is smaller than the cutoff angle
    if (diff_angle <= (p.cross_severing_start_angle) * PI / 180.)
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
    cut_seg_number = random_gen.rand_int(cut_is->second.occupancy - 1); // nummers 0 - occupancy -1 .
    if (cut_seg_number < cut_is->second.occupancy / 2)                // search from the beginning.
                                                    // In principle, all searches could start from the beginning. Also
                                                    // starting from the end should save (on average) half of the time
    {
        test_seg = cut_is->second.mirror->second.other_trajectory->segments.begin();
        final_test_seg = cut_is->second.mirror->second.other_trajectory->segments.end();
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
        test_seg = cut_is->second.mirror->second.other_trajectory->segments.end();
        final_test_seg = cut_is->second.mirror->second.other_trajectory->segments.begin();
        while (test_seg != final_test_seg)
        {
            test_seg--;
            if ((**test_seg).crosses_intersection(cut_is))
            {
                cut_seg_number++;
                if (cut_seg_number == cut_is->second.occupancy)
                {
                    cut_seg = *test_seg;
                    break;
                }
            }
        }
    }
    if (cut_seg != nullptr)
    {
        cut_mt = cut_seg->mt;
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
    if ((p.restricted_pool != 0) && (p.v_tm / (p.v_plus * (1. - total_length / (p.pool_density * geometry->area))) > 0.9999))
    {
        return;
    }

    // no nucleation at forbidden zone
    TrajectoryVector tv;

    // get random nucleation point (vector)
    SurfaceVector sv = geometry->random_surface_vector();

    // forbiden zone and parameter modified : no microtubule created
    if (sv.region->face_tag == -1)
    {
        return;
    }

    double override_angle(0.0);

    switch (p.nucleation_type)
    {
            // select isotropic nucleation
        case NucleationType::ISOTROPIC:
            override_angle = sv.angle;
            break;
        default:
            break;
    }

    sv.angle = override_angle;
    tv = geometry->create_trajectory(sv);
    growing_mts.create(this, tv);

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

    if ((growing_mts.size() == 0) || (p.restricted_pool == 0))
    {
        double const_rate(0.0);
        double slope_rate(0.0);

        const_rate = p.k_nuc * geometry->area + p.k_res * shrinking_mts.size() + p.k_sev * total_length;

        for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
        {
            const_rate += p.k_cat * p.region_kcat_multiplier[d_tag] * growing_tips_reg[d_tag];
        }

#ifdef CROSS_SEV
        const_rate += p.k_cross * occupied_intersection_list.size();
#endif

        slope_rate = p.k_sev * ((p.v_plus - p.v_tm) * growing_mts.size() + (p.v_min - p.v_tm) * shrinking_mts.size());

        // generate a random number larger than 0 with an exponential distribution
        log_random = -log(random_gen.rand_dbl_exc());

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
        double inv_l_max = 1.0 / (p.pool_density * geometry->area);
        double alpha = (p.v_plus - p.v_tm) * growing_mts.size() + (p.v_min - p.v_tm) * shrinking_mts.size();
        double beta = p.v_plus * growing_mts.size() * inv_l_max;
        double l_inf = alpha / beta;

        double base_rate = p.k_nuc * geometry->area + p.k_res * shrinking_mts.size();
        for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
        {
            base_rate += p.k_cat * p.region_kcat_multiplier[d_tag] * growing_tips_reg[d_tag];
        }

#ifdef CROSS_SEV
        base_rate += p.k_cross * occupied_intersection_list.size();
#endif

        double x = p.k_sev * (total_length - l_inf) / (base_rate + p.k_sev * l_inf);
        double y = -beta * log(random_gen.rand_dbl_exc()) / (base_rate + p.k_sev * l_inf);

        next_event_interval = (-x + y + lambert_w(x * exp(x - y))) / beta;

        double total_rate = base_rate + p.k_sev * (l_inf + (total_length - l_inf) * exp(-beta * next_event_interval));
        inv_total_rate = 1.0 / total_rate;
    }

#ifdef DBG_ACID_TEST
    if ((nextEventInterval < 0) || (nextEventInterval != nextEventInterval))
    {
        cerr << "Invalid deterministic event time calculated. Exiting\n";
        exit(-1);
    }
#endif

    next_stochastic_event_time = system_time + next_event_interval;

    // even if  invTotalRate is negative! due to the timing, the event will never be executed.
    type_selector = random_gen();

    double type_select_reg_source(0.0);
    for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
    {
        type_select_reg_source += p.k_cat * p.region_kcat_multiplier[d_tag] * growing_tips_reg[d_tag];
    }

    if (type_selector < p.k_nuc * geometry->area * inv_total_rate)
    {
        next_stochastic_event_type = StochasticEventType::NUCLEATION;
    }

    // select catastrophe event
    else if (type_selector < (p.k_nuc * geometry->area + type_select_reg_source) * inv_total_rate)
    {
        next_stochastic_event_type = StochasticEventType::CATASTROPHE;
    }

    // select rescue event
    else if (type_selector
             < (p.k_nuc * geometry->area + type_select_reg_source + p.k_res * shrinking_mts.size()) * inv_total_rate)
    {
        next_stochastic_event_type = StochasticEventType::RESCUE;
    }

// select severing event
#ifdef CROSS_SEV
    else if (type_selector < (p.k_nuc * geometry->area + type_select_reg_source + p.k_res * shrinking_mts.size()
                             + p.k_cross * occupied_intersection_list.size())
                            * inv_total_rate)
    {
        next_stochastic_event_type = StochasticEventType::SEVERING_AT_CROSS;
    }
#endif

    // select katanin induced event
    else
    {
        next_stochastic_event_type = StochasticEventType::KATANIN;
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

    switch (event.global_type)
    {

        case GlobalEventType::MEASURE:
            break;
        case GlobalEventType::STATUS:
            if (!integrity_check())
            {
                cerr << "Integrity check failed. Exiting.\n";
                exit(-1);
            }

            perform_measurement();
            next_status_event_time = system_time + p.measurement_interval;
            time_queue.push_global(next_status_event_time, GlobalEventType::STATUS);
            break;

        case GlobalEventType::SNAPSHOT:
            // declaring variable within case: so the local blocking
            {
                OrderParameters order;

                geometry->get_order_parameters(order);

                // write a movie frame
                movie_file << order.rdirector[0] << "," << order.rdirector[1] << "," << order.rdirector[2] << "#";
                movie_file << geometry->object_cm(0, 0) << "," << geometry->object_cm(1, 0) << ","
                          << geometry->object_cm(2, 0) << "#";
                movie_file << geometry->object_pa(0, 0) << "," << geometry->object_pa(1, 0) << ","
                          << geometry->object_pa(2, 0) << "*";

                geometry->output_snapshot(movie_file);

                // write a local order parameter heatmap frame
                geometry->output_order_heat_map(heat_map_file, order.local_order, order.sv);

                next_snapshot_event_time = system_time + p.movie_frame_interval;
                time_queue.push_global(next_snapshot_event_time, GlobalEventType::SNAPSHOT);
            }
            break;

        case GlobalEventType::STOP:
        {
            OrderParameters order;
            geometry->get_order_parameters(order);

            // division plane: polygon edge-list & area
            Eigen::Vector3d normal_ppb(order.rdirector[0], order.rdirector[1], order.rdirector[2]);
            geometry->area_ppb = intersecting_polygon(geometry->ppb,
                                                    geometry->global_vertex,
                                                    geometry->regions,
                                                    geometry->ppb_edge_list,
                                                    normal_ppb,
                                                    geometry->nucleous_position,
                                                    p.geometry);

            if (p.show_output == 1)
            {
                cout << "Division plane area: " << geometry->area_ppb << endl;
            }

            stop_signal = true;
        }
        break;

        case GlobalEventType::PARAMETER_CHANGE:
        {
            cout << "Changing parameters. Reading from file [" << p.new_parameter_file << "]\n";
            flush_and_reload(false);
            string old_dir = p.output_dir;
            if (p.reinitialize(p.new_parameter_file.c_str()))
            {
                if (p.output_dir != old_dir)
                {
                    cout << "Switching directories.\n";
                    write_measurements_to_file(0);
                    close_files();
                    initialize_output();
                }
                p.write_to_file();

                if (p.movie_enabled)
                {
                    if (next_snapshot_event_time < ZERO_CUTOFF)
                    {
                        time_queue.push_global(p.movie_frame_interval, GlobalEventType::SNAPSHOT);
                        next_snapshot_event_time = p.movie_frame_interval;
                    }
                }
                else
                {
                    next_snapshot_event_time = 0;
                }

                if (p.new_parameter_read_interval > ZERO_CUTOFF)
                {
                    time_queue.push_global(p.new_parameter_read_interval, GlobalEventType::PARAMETER_CHANGE);
                    next_parameter_event_time = p.new_parameter_read_interval;
                }
            }
            else
            {
                cerr << "Failed to read parameters from file. Exiting.\n";
                emergency_break();
                exit(-1);
            }

            time_queue.push_global(p.new_parameter_read_interval, GlobalEventType::PARAMETER_CHANGE);
            next_parameter_event_time = p.new_parameter_read_interval;
            break;
        }

        case GlobalEventType::SIGNAL_PPB:
            cout << " Starting formation of PPB ...\n";
            p.write_to_file();
            p.pp_bk_nuc_fraction = 0.8;
            cout << p.pp_bk_nuc_fraction << " : " << p.ppb << endl;
            OrderParameters order;
            geometry->get_order_parameters(order);
            Eigen::Vector3d normal_ppb(order.rdirector[0], order.rdirector[1], order.rdirector[2]);
            geometry->area_ppb = intersecting_polygon(geometry->ppb,
                                                    geometry->global_vertex,
                                                    geometry->regions,
                                                    geometry->ppb_edge_list,
                                                    normal_ppb,
                                                    geometry->nucleous_position,
                                                    p.geometry);

            for (int d_tag = 0; d_tag <= p.face_number; d_tag++)
            {
                growing_tips_reg[d_tag] = 0;
            }

            for (size_t eno = 0; eno < geometry->regions.size(); eno++)
            {
                if (geometry->regions[eno]->poly_intersect_mark == 1)
                {
                    geometry->regions[eno]->face_tag = 1;

                    for (size_t tip = 0; tip < geometry->regions[eno]->growing_plus_tip_list.size(); tip++)
                    {
                        growing_tips_reg[1] += 1;
                    }

                    geometry->patch_area[1] += geometry->regions[eno]->area;
                    geometry->regions_index[1].push_back(eno);
                }

                else
                {
                    geometry->regions[eno]->face_tag = 2;

                    for (size_t tip = 0; tip < geometry->regions[eno]->growing_plus_tip_list.size(); tip++)
                    {
                        growing_tips_reg[2] += 1;
                    }

                    geometry->patch_area[2] += geometry->regions[eno]->area;
                    geometry->regions_index[2].push_back(eno);
                }
            }
            time_queue.push_global(VERY_LARGE, GlobalEventType::SIGNAL_PPB);
            next_ppb_event_time = VERY_LARGE;

            break;
    }

    return;
}

EventDescriptorIndex System::register_event_descriptor(EventDescriptor* ei)
{
    // maybe use a pool of numbers instead of an increasing number....this could wrap around
    EventDescriptorIndex index = event_descriptor_id++;
    event_descriptor_map[index] = ei;

    if (event_descriptor_id == numeric_limits<EventDescriptorIndex>::max() - 10)
    {
        time_queue.push_global(system_time, GlobalEventType::STOP);
    }

    return index;
}

void System::unregister_event_descriptor(EventDescriptorIndex ei)
{
    event_descriptor_map.erase(event_descriptor_map.find(ei));
    return;
}

EventDescriptor* System::get_event_descriptor(EventDescriptorIndex idx)
{
    map<EventDescriptorIndex, EventDescriptor*>::iterator itr = event_descriptor_map.find(idx);
    if (itr == event_descriptor_map.end())
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
    occupied_intersection_list.remove_element(is.occupied_list_ptr->index);
    is.occupied_list_ptr = nullptr;
    is.mirror->second.occupied_list_ptr = nullptr;
    return;
}

void System::add_occupied_intersection(IntersectionItr is)
{
    is->second.occupied_list_ptr = occupied_intersection_list.create(is);
    is->second.mirror->second.occupied_list_ptr = is->second.occupied_list_ptr;
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
            binomial_table[i][j] = factorial[i] / (factorial[j] * factorial[i - j]);
        }
    }
    return;
}
