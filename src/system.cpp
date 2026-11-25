#include "mesh.hpp"
#include "region.hpp"
#include "system.hpp"
#include "mt_tip.hpp"
#include "parameters.hpp"

void System::collision_probabilities(double angle, double& pcat, double& pzip)
{
    switch (p.interactionType)
    {
        case INT_ZIP_FIRST:

            // angle < magic angle
            if (angle < PI * p.magicAngle / 180.0)
            {
                pcat = 0.0;

                // zippering enabled
                if (p.zipperingEnabled)
                {
                    pzip = p.zipFraction;
                }
                else
                {
                    pzip = 0.0;
                }
                break;
            }

            pzip = 0.0;

            // angle > p.catStartAngle and induced catastrophes enabled
            if (p.catastrophesEnabled && (angle > p.catStartAngle / 180.0 * PI))
            {
                // proportional catastrophe enabled
                if (p.proportionalCatastrophes)
                {
                    pcat = (angle - p.catStartAngle / 180.0 * PI) / (PI / 2. - p.catStartAngle / 180.0 * PI)
                           * p.inducedCatastropheFraction;
                }

                // proportional catastrophe disabled
                else
                {
                    pcat = p.inducedCatastropheFraction;
                }
            }
            else
            {
                pcat = 0.0;
            }

            break;

        case INT_CAT_FIRST:

            // angle > p.catStartAngle and induced catastrophes enabled
            if (p.catastrophesEnabled && (angle > p.catStartAngle / 180.0 * PI))
            {
                // proportional catastrophe enabled
                if (p.proportionalCatastrophes)
                {
                    pcat = (angle - p.catStartAngle / 180.0 * PI) / (PI / 2. - p.catStartAngle / 180.0 * PI)
                           * p.inducedCatastropheFraction;
                }

                // proportional catastrophe disabled
                else
                {
                    pcat = p.inducedCatastropheFraction;
                }
            }
            else
            {
                pcat = 0.0;
            }

            // angle < magic angle and zippering enabled
            if (p.zipperingEnabled && (angle < p.magicAngle / 180.0 * PI))
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
    double totalSum = 0;

    int i, j;

    for (i = 0; i <= npar; i++)
    {
        for (j = max(0, i - ncoll); j <= max(0, i - 1); j++)
        {
            totalSum += (npar + 1 - i) * binomialTable[ncoll][i - j] * binomialTable[max(0, i - 1)][j]
                        * pow(x_single, 2 * j + ncoll - i) * pow(z_single, 2 * (i - j));
        }
    }
    return totalSum / (npar + 1);
}

double System::multi_pzip(int n_par, int n_coll, double x_single, double z_single)
{
    double totalSum = 0;

    int i, j, k;

    for (i = 0; i <= n_par; i++)
    {
        for (j = 0; j <= n_coll - 1; j++)
        {
            for (k = max(0, i - j); k <= i; k++)
            {
                totalSum += binomialTable[j][i - k] * binomialTable[i][k] * pow(x_single, 2 * k + j - i)
                            * pow(z_single, 2 * (i - k) + 1);
            }
        }
    }
    return totalSum / (n_par + 1);
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

    double Pcat{ 0.0 };
    double Pzip{ 0.0 };
    double Psum{ 0.0 };

    // get a random number
    double rand = randomGen.rand_dbl_exc();

    // get the available Pcat and Pzip
    collision_probabilities(angle, Pcat, Pzip);

    switch (p.bundleType)
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
                if (rand < Pcat)
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
            if (rand < Pcat)
            {
                return CT_INDUCED_CATASTROPHE;
            }

            // zippering selected
            else if (rand < Pcat + Pzip)
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
            Pcat *= n_cross;
            Pzip *= n_cross;
            Psum = Pcat + Pzip;

            // now normalize the probabilities
            if (Psum > 1)
            {
                Pcat /= Psum;
                Pzip /= Pzip;
            }

            if (rand < Pcat)
            {
                return CT_INDUCED_CATASTROPHE;
            }

            else if (rand < Pcat + Pzip)
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

            double Pcross = 1.0 - Pcat - Pzip;
            double multiPx = multi_pcross(min(n_parallel, MAXBINOM - 1), min(n_cross, MAXBINOM - 1), Pcross, Pzip);
            double multiPz = multi_pzip(min(n_parallel, MAXBINOM - 1), min(n_cross, MAXBINOM - 1), Pcross, Pzip);

            if (rand < multiPx)
            {
                return CT_CROSSOVER;
            }

            else if (rand < multiPx + multiPz)
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
    geometry(nullptr),
    p(this),
    systemTime(0),
    systemTimeOffset(0),
    currentTimeTag(0),
    totalLength(0),
    memUsage(0),
    stopSignal(false),
    wallClockStartTime(time(nullptr)),
    totalSEventCount(0),
    totalValidDEventCount(0),
    totalInvalidDEventCount(0),
    totalLengthSeveringCount(0),
    totalIntersectionSeveringCount(0),
    boundaryCrossingCount(0),
    totalZipperCount(0),
    totalCrossoverCount(0),
    totalInducedCatastropheCount(0),
    countSegments(0),
    countTrajectories(0),
    countIntersections(0),
    timeQueue(this, &System::identity, &System::identity),
    vPlusQueue(this, &System::v_plus_to_time, &System::time_to_v_plus),
    nextStatusEventTime(0),
    nextSnapshotEventTime(0),
    EventDescriptorID(0),
    eventID(0),
    dataSaved(false)
{
    /// read parameters
    p.initialize(par_file);

    // initialize number of growing tips in different faces
    for (int dTag = 0; dTag <= p.faceNumber; dTag++)
    {
        growingTipsReg.push_back(0);
    }

    // load geometry
    geometry = new TriMeshGeometry(this);

    /// initialize the random seed with the desired one
    randomGen.seed(p.seed);

    /// used in multiPcross & multiPzip events
    make_binomial_table();

    /// perform further initialization
    initialize_output();
    p.write_to_file();

    /// start creating MT events
    determine_stochastic_event();
    timeQueue.push_global(p.measurementInterval, STATUS);
    nextStatusEventTime = p.measurementInterval;

    /// if snapShot events are requested
    if (p.movieEnabled)
    {
        timeQueue.push_global(p.movieFrameInterval, SNAPSHOT);
        nextSnapshotEventTime = p.movieFrameInterval;
    }

    /// if parameter file is updated on run
    if (p.PPB)
    {
        timeQueue.push_global(p.PPBformingTime, SIGNAL_PPB);
        nextPPBEventTime = p.PPBformingTime;
    }

    /// if parameter file is updated on run
    if (p.newParameterReadInterval > ZERO_CUTOFF)
    {
        timeQueue.push_global(p.newParameterReadInterval, PARAMETER_CHANGE);
        nextParameterEventTime = p.newParameterReadInterval;
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
    int wallTimePollingCounter(0);
    int memoryPollingCounter(0);
    timeQueue.push_global(sim_time - systemTimeOffset, STOP);

    if (p.showOutput == 1)
    {
        cout << "===============================================" << endl;
        cout << " Time stamped simulation output " << endl;
        cout << "===============================================" << endl;
    }

    while (!stopSignal)
    {
        if (++wallTimePollingCounter == CLOCK_POLLING_INTERVAL)
        {
            if (difftime(time(nullptr), wallClockStartTime) >= p.wallClockLimit)
            {
                // time's up; send a stop signal
                cerr << "Wall clock limit (" << p.wallClockLimit << " seconds) reached. Stopping.\n";
                timeQueue.push_global(systemTime, STOP);

                // write final snapshot when running out of time.
                if (p.movieEnabled)
                {
                    geometry->output_snapshot(movieFile);
                }

                // emergency exit the program
                emergency_break();
                exit(-1);
            }
            wallTimePollingCounter = 0;
        }

        if (++memoryPollingCounter == MEMORY_POLLING_INTERVAL)
        {
            if (estimate_memory_footprint() > 1024 * 1024 * p.memoryLimit)
            {
                cerr << "Memory limit (approx. " << p.memoryLimit << " MB) reached. Stopping.\n";
                timeQueue.push_global(systemTime, STOP);

                // emergency exit the program
                emergency_break();
                exit(-1);
            }
            memoryPollingCounter = 0;
        }

        if (eventID > QUEUE_FLUSH_INTERVAL)
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
    total += EventDescriptorMap.size() * (sizeof(pair<EventDescriptorIndex, EventDescriptor*>) + 3 * sizeof(void*));
    total += countSegments * (sizeof(Segment) + sizeof(Segment*));
    total += countTrajectories * (sizeof(Trajectory));
    total += countIntersections * (2 * (sizeof(pair<double, Intersection>) + 3 * sizeof(void*)));
    total += (timeQueue.queue.size() + vPlusQueue.queue.size()) * (sizeof(DeterministicEvent) + 2 * sizeof(void*));
    total += measurementHistory.size() * sizeof(Measurement);

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
    double timeOffset(systemTime);

    // update all microtubule lengths
    update_all();

    // recalculate total lengths to avoid drift (in DBG mode, also check for large deviations)
    double regionLength(0.);
    double systemLength(0);
    Trajectory* tptr(nullptr);

    for (size_t ridx = 0; ridx < geometry->regions.size(); ridx++)
    {
        regionLength = 0;
        tptr = geometry->regions[ridx]->trajectories.first();

        while (tptr != nullptr)
        {
            regionLength += tptr->segment_length();
            tptr = tptr->next();
        }

        if (abs(geometry->regions[ridx]->totalLength - regionLength) > ZERO_CUTOFF * max(100., regionLength))
        {
            cerr << "Unacceptable drift in region length (drift = "
                 << abs(geometry->regions[ridx]->totalLength - regionLength)
                 << "). Exiting. [try lowering flush interval]\n";
            emergency_break();
            exit(-1);
        }
        geometry->regions[ridx]->totalLength = regionLength;
        systemLength += regionLength;
    }

    if (abs(totalLength - systemLength) > ZERO_CUTOFF * max(1., totalLength))
    {
        cerr << "Unacceptable drift in system length (drift = " << abs(totalLength - systemLength)
             << "). Exiting. [try lowering flush interval]\n";
        emergency_break();
        exit(-1);
    }

    totalLength = systemLength;

    // reset master timer
    systemTimeOffset += timeOffset;

    systemTime = 0;
    timeQueue.flush();
    vPlusQueue.flush();
    EventDescriptorMap.clear();
    eventID = 0;
    EventDescriptorID = 0;

    // now, reload all events (this is a PROBLEM for global events that have not been saved. They need to be re-inserted
    // manually)
    determine_stochastic_event();
    timeQueue.push_global(p.stopTime - systemTimeOffset, STOP);
    nextStatusEventTime -= timeOffset;
    timeQueue.push_global(nextStatusEventTime, STATUS);

    if (refresh_parameter_event)
    {
        nextParameterEventTime -= timeOffset;
        timeQueue.push_global(nextParameterEventTime, PARAMETER_CHANGE);
    }

    if (p.movieEnabled)
    {
        nextSnapshotEventTime -= timeOffset;
        timeQueue.push_global(nextSnapshotEventTime, SNAPSHOT);
    }

    if (p.PPB)
    {
        nextPPBEventTime -= timeOffset;
        timeQueue.push_global(nextPPBEventTime, SIGNAL_PPB);
    }

    Microtubule* mt(nullptr);
    Segment* seg(nullptr);
    mt = growing_mts.first();

    while (mt != nullptr)
    {
        // this can go into a change time function
        mt->nucleationTime -= timeOffset;
        seg = mt->segments.first();

        while (seg != nullptr)
        {
            seg->nucleationTime -= timeOffset;
            seg = seg->next();
        }

        mt->plus.event.index = register_event_descriptor(&(mt->plus.event));
        mt->minus.event.index = register_event_descriptor(&(mt->minus.event));
        mt->disappearEvent.index = register_event_descriptor(&(mt->disappearEvent));
        mt->plus.determine_event();
        mt->minus.determine_event();
        mt->set_disappear_event();
        mt = mt->next();
    }

    mt = shrinking_mts.first();

    while (mt != nullptr)
    {
        mt->nucleationTime -= timeOffset;
        seg = mt->segments.first();

        {
            seg->nucleationTime -= timeOffset;
            seg = seg->next();
        }

        mt->plus.event.index = register_event_descriptor(&(mt->plus.event));
        mt->minus.event.index = register_event_descriptor(&(mt->minus.event));
        mt->disappearEvent.index = register_event_descriptor(&(mt->disappearEvent));
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
    const double expMin1 = exp(-1.0);
    const double zeroCutoff = pow(2.0, -DBL_MANT_DIG / 6.0);
    const double etaCutoff = 34.0 * pow(2.0, 2.0 * (1.0 - DBL_MANT_DIG) / 7.0);

    // LambertW function (Barry et al., ACM Trans. on Math. Software, Vol. 21, No. 2, June 1995, pp. 161--171)
    if (x < -expMin1)
    {
        return -VERY_LARGE;
    }

    if (abs(x) < zeroCutoff)
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
        double D(0.0);
        if (sqrt_eta >= etaCutoff)
        {
            double n2 = 4.612634277343749 * sqrt(sqrt(sqrt_eta + 1.09556884765625));
            double n1 = (4.0 - 3.0 * sqrt(2.0) + (2.0 * sqrt(2.0) - 3.0) * n2) / (sqrt(2.0) - 2.0);
            D = n1 * sqrt_eta / (n2 + sqrt_eta);

            out = -1.0 + sqrt_eta / (1.0 + sqrt_eta / (3.0 + D));

            double z_n = log(x / out) - out;
            double temp = 2.0 * (1.0 + out) * (1.0 + out + z_n * 2.0 / 3.0);
            double e_n = (z_n / (1.0 + out)) * (temp - z_n) / (temp - 2.0 * z_n);
            out = out * (1.0 + e_n);
        }
        else
        {
            D = sqrt_eta
                / (8.0 / 3.0 + sqrt_eta / (135.0 / 83.0 + sqrt_eta / (166.0 / 39.0 + sqrt_eta * 3167.0 / 3549.0)));
            out = -1.0 + sqrt_eta / (1.0 + sqrt_eta / (3.0 + D));
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
    if (p.restrictedPool == 0)
    {
        return dist;
    }

    // for finite tubulin pool
    double Lmax = (p.poolDensity * geometry->area);
    double invLMax = 1.0 / (p.poolDensity * geometry->area);

    double alpha = (p.vPlus - p.vTM) * growing_mts.size() + (p.vMin - p.vTM) * shrinking_mts.size();

    if (growing_mts.size() == 0)
    {
        if (shrinking_mts.size() == 0)
        {
            return dist / (1.0 - totalLength * invLMax);
        }
        else
        {
            return ((Lmax - totalLength) - sqrt(pow(Lmax - totalLength, 2) - 2. * alpha * Lmax * dist)) / alpha;
        }
    }

    double beta = p.vPlus * growing_mts.size() * invLMax;

    if ((shrinking_mts.size() == 0) && ((p.vTM / p.vPlus) < ZERO_CUTOFF))
    {
        return -log(1.0 - growing_mts.size() * p.vPlus * dist / (Lmax - totalLength)) / beta;
    }

    double pv = (alpha / beta - totalLength) / (Lmax - alpha / beta);
    double qv = growing_mts.size() * p.vPlus * dist / (Lmax - alpha / beta);

    // verify whether this is an optimal trade-off
    if (pv - qv < 250)
    {
        return (qv - pv + lambert_w(pv * exp(pv - qv))) / beta;
    }

    // if not, do asymptotic approximation
    double L1 = pv - qv + log(pv);
    double L2 = log(L1);
    return (1.0 / beta)
           * (qv - pv + L1 - L2 + L2 / L1 + L2 * (L2 - 2) / (2 * pow(L1, 2))
              + L2 * (6 - 9 * L2 + 2 * pow(L2, 2)) / (6 * pow(L1, 3))
              + L2 * (-12 + 36 * L2 - 22 * pow(L2, 2) + 3 * pow(L2, 3)) / (12 * pow(L1, 4)));
}

double System::time_to_v_plus(double time)
{
    // for infinite pool (implimented as identity function)
    if (p.restrictedPool == 0)
    {
        return time;
    }

    // for finite tubulin pool
    double invLMax = 1.0 / (p.poolDensity * geometry->area);
    double alpha = (p.vPlus - p.vTM) * growing_mts.size() + (p.vMin - p.vTM) * shrinking_mts.size();

    if (growing_mts.size() == 0)
    {
        return ((1.0 - totalLength * invLMax) - 0.5 * alpha * invLMax * time) * time;
    }

    double beta = p.vPlus * growing_mts.size() * invLMax;
    double lambda_inf = invLMax * alpha / beta;
    double lambda_0 = totalLength * invLMax;

    return (1.0 - lambda_inf) * time - (lambda_0 - lambda_inf) * (1.0 - exp(-beta * time)) / beta;
}

void System::advance_time(double new_time)
{
#ifdef DBG_SYSTEM
    cout << "DBG/SYSTEM: advancing time by " << newTime - systemTime << " seconds. New system time is " << newTime
         << ".\n";
#endif

    // make a jump in time step (taking advantage of using event-driven algorithm)
    timeQueue.advance_time(new_time);
    vPlusQueue.advance_time(new_time);
    systemTime = new_time;
    totalLength += timeQueue.progression(currentTimeTag)
                   * (p.vMin * shrinking_mts.size() - p.vTM * (growing_mts.size() + shrinking_mts.size()))
                   + vPlusQueue.progression(currentTimeTag) * (p.vPlus * growing_mts.size());

    if (++currentTimeTag == POSITION_CACHE_SIZE)
    {
        currentTimeTag = 0;
        update_all(true);
    }

    timeQueue.store_time(currentTimeTag);
    vPlusQueue.store_time(currentTimeTag);

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
    double tqValue = numeric_limits<double>::max();
    double dqValue = numeric_limits<double>::max();

    // maybe clean up a little
    if (!timeQueue.empty())
    {
        tqValue = timeQueue.first_event_time();
    }

    if (!vPlusQueue.empty())
    {
        dqValue = vPlusQueue.first_event_time();
    }

    if ((nextStochasticEventTime < tqValue) && (nextStochasticEventTime < dqValue))
    {
        // execute stochastic event
        totalSEventCount++;
#ifdef DBG_EVENT
        cout << "DBG/Event: executing stochastic event [time=" << nextStochasticEventTime
             << ", type=" << nextStochasticEventType << "]\n";
#endif

        advance_time(nextStochasticEventTime);

        // execute next stochastic event
        switch (nextStochasticEventType)
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
        double nextEventTime(0.0);

        // execute deterministic event
        DeterministicEvent event;
        if (tqValue < dqValue)
        {
            event = timeQueue.pop();
            nextEventTime = tqValue;
        }
        else
        {
            event = vPlusQueue.pop();
            nextEventTime = dqValue;
        }

        if (event.infoIdx != -1)
        {
            map<EventDescriptorIndex, EventDescriptor*>::iterator eventItr = EventDescriptorMap.find(event.infoIdx);

            // whether the tip still exists and the tag still matches (valid event)
            if ((eventItr != EventDescriptorMap.end()) && (eventItr->second->tag == event.tag))
            {
#ifdef DBG_EVENT
                cout << "DBG/Event: executing deterministic event [time=" << event.eventTimeDist
                     << ", type=" << eventItr->second->type << ", tag=" << event.tag << "]\n";
#endif

                totalValidDEventCount++;
                advance_time(nextEventTime);
                eventItr->second->mt->handle_event(eventItr->second);
                determine_stochastic_event();
            }

            // invalid event
            else
            {
                totalInvalidDEventCount++;
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

            advance_time(event.eventTimeDist);
            handle_global_event(event);
        }
    }

    return;
}

void System::handle_catastrophe_event()
{

    size_t tipNumber(0), totalRelevantTips(0), ridx(0);
    int domainType(0);
    RegionMTTipTag tipTag;

    double norm(0.0);
    for (int dTag = 0; dTag <= p.faceNumber; dTag++)
    {
        norm += growingTipsReg[dTag] * p.RegionKcatMultiplier[dTag];
    }

    double randGen = randomGen.rand();

    double normFraction(0.0);
    for (int dTag = 0; dTag <= p.faceNumber; dTag++)
    {
        normFraction += growingTipsReg[dTag] * p.RegionKcatMultiplier[dTag];
        if (randGen < (double)(normFraction) / norm)
        {
            domainType = dTag;
            totalRelevantTips = growingTipsReg[dTag];
            break;
        }
    }

    tipNumber = randomGen.rand_int(totalRelevantTips - 1);

    // determine the region this tip is located in
    // depending on whether it's a high or low number, start scanning from the beginning or end
    if (tipNumber < totalRelevantTips / 2)
    {
        while (true)
        {
            if (geometry->regions[ridx]->faceTag == domainType)
            {
                if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size())
                {
                    break;
                }

                tipNumber -= geometry->regions[ridx]->growingPlusTipList.size();
            }
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tipNumber = totalRelevantTips - 1 - tipNumber;
        ridx = geometry->regions.size() - 1;
        while (true)
        {
            if (geometry->regions[ridx]->faceTag == domainType)
            {
                if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size())
                {
                    break;
                }
                tipNumber -= geometry->regions[ridx]->growingPlusTipList.size();
            }
            ridx--;
        }
        // flip the direction of the result for the next step
        tipNumber = geometry->regions[ridx]->growingPlusTipList.size() - 1 - tipNumber;
    }

    // Now, select the tip with # 'tipNumber' from the current region.
    // Again, if the number is high within the number of tips in the region, start from the end
    if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size() / 2)
    {
        tipTag = geometry->regions[ridx]->growingPlusTipList.begin();
        while (tipNumber > 0)
        {
            tipTag++;
            tipNumber--;
        }
    }
    else
    {
        tipTag = --(geometry->regions[ridx]->growingPlusTipList.end());
        while (tipNumber < geometry->regions[ridx]->growingPlusTipList.size() - 1)
        {
            tipTag--;
            tipNumber++;
        }
    }

    (*tipTag)->mt->catastrophe();

    return;
}

void System::handle_rescue_event()
{
    // TODO
    // Fix the logic below and the types of these variables
    // to eliminate the need for static_cast.
    size_t ridx(0);
    unsigned long int tipNumber(0);

    RegionMTTipTag tipTag;

    // select the number of the tip to be rescued
    tipNumber = randomGen.rand_int(shrinking_mts.size() - 1);

    // determine the region this tip is located : binary search algorithm
    if (tipNumber < static_cast<unsigned long long>(shrinking_mts.size() / 2))
    {
        while (tipNumber >= geometry->regions[ridx]->shrinkingPlusTipList.size())
        {
            tipNumber -= geometry->regions[ridx]->shrinkingPlusTipList.size();
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tipNumber = shrinking_mts.size() - 1 - tipNumber;
        ridx = geometry->regions.size() - 1;
        while (tipNumber >= geometry->regions[ridx]->shrinkingPlusTipList.size())
        {
            tipNumber -= geometry->regions[ridx]->shrinkingPlusTipList.size();
            ridx--;
        }
        // flip the direction of the result for the next step
        tipNumber = geometry->regions[ridx]->shrinkingPlusTipList.size() - 1 - tipNumber;
    }

    // determine which tip in the selected region : binary search algorithm
    if (tipNumber < geometry->regions[ridx]->shrinkingPlusTipList.size() / 2)
    {
        tipTag = geometry->regions[ridx]->shrinkingPlusTipList.begin();
        while (tipNumber > 0)
        {
            tipTag++;
            tipNumber--;
        }
    }
    else
    {
        tipTag = --(geometry->regions[ridx]->shrinkingPlusTipList.end());
        while (tipNumber < geometry->regions[ridx]->shrinkingPlusTipList.size() - 1)
        {
            tipTag--;
            tipNumber++;
        }
    }

    // rescue the selected tip
    (*tipTag)->mt->rescue();

    return;
}

void System::handle_severing_event()
{
    double cutLength(0.0);
    Segment* cutSeg(nullptr);
    random_position_on_microtubule(cutLength, cutSeg);

    // sever the MT now !
    if (cutSeg != nullptr)
    {
        cutSeg->mt->sever(cutSeg, cutLength);
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
    double cutLength = randomGen.rand_dbl_exc(totalLength);
    // first, select relevant region. Go from both sides to reduce time.
    if (cutLength < 0.5 * totalLength)
    {
        ridx = 0;
        while (true)
        {
            geometry->regions[ridx]->update_region_length();
            if (cutLength < geometry->regions[ridx]->totalLength)
            {
                break;
            }
            cutLength -= geometry->regions[ridx]->totalLength;
            ridx++;
            if (ridx == static_cast<long long>(geometry->regions.size()))
            {
                // apparently, cutLength was *just* too long...
                cutLength -= ZERO_CUTOFF;
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
        cutLength = totalLength - cutLength;
        while (true)
        {
            geometry->regions[ridx]->update_region_length();
            if (cutLength < geometry->regions[ridx]->totalLength)
            {
                break;
            }
            cutLength -= geometry->regions[ridx]->totalLength;
            ridx--;
            if (ridx == -1)
            {
                // apparently, cutLength was *just* too long...
                cutLength -= ZERO_CUTOFF;

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
        cutLength = geometry->regions[ridx]->totalLength - cutLength;
    }

    // Now that the area has been selected, cycle over trajectories and their segments
    Trajectory* trptr(nullptr);
    list<Segment*>::iterator seg;
    double temp(0.0);
    bool stopFlag(false);
    // Again, let the search direction depend on the expected position
    if (cutLength < 0.5 * (geometry->regions[ridx]->totalLength))
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
                if (cutLength < (temp = (*seg)->length()))
                {
                    stopFlag = true;
                    break;
                }
                else
                {
                    // if not, continue
                    cutLength -= temp;
                    seg++;
                }
            }

            if (stopFlag)
            {
                break;
            }

            // if we ran out of trajectories, something is wrong. This could happen due to numerical inaccuracies.
            if (trptr->next() == nullptr)
            {
                cerr << "Rare event: cutting on the edge\n";
                seg--;
                cutLength -= ZERO_CUTOFF;
                break;
            }
            trptr = trptr->next();
        }
    }
    else
    {
        // flip the length and trajectory search orders
        cutLength = geometry->regions[ridx]->totalLength - cutLength;
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
                if (cutLength < (temp = (**seg).length()))
                {
                    stopFlag = true;
                    break;
                }
                else
                {
                    cutLength -= temp;
                }
            } while (seg != trptr->segments.begin());

            if (stopFlag)
            {
                break;
            }

            if (trptr->previous() == nullptr)
            {
                cerr << "Rare event: cutting on the edge\n";
                seg++;
                cutLength -= ZERO_CUTOFF;
                break;
            }
            trptr = trptr->previous();
        }
        cutLength = (*seg)->length() - cutLength;
    }

    random_seg = *seg;
    random_pos = cutLength;
    return;
}

#ifdef CROSS_SEV

// TODO: verify this function
void System::handle_severing_at_cross_event()
{
    Microtubule* cutMT(nullptr);
    IntersectionItr cutIS;
    OccupiedIntersection* cutOccIS(nullptr);
    Segment* cutSeg(nullptr);
    int cutSegNumber(0);
    list<Segment*>::iterator testSeg;
    list<Segment*>::iterator finalTestSeg;

    cutOccIS = OccupiedIntersectionList.element_address(randomGen.rand_int(OccupiedIntersectionList.size() - 1));
    if (p.crossSeveringTop || (randomGen.rand_int(1) == 0))
    {
        // if MT (bundle) on top should be cut, select that one, or (if random) pick it with 50% chance.
        // Note that no effort is made to ensure proportional selection for bundles of different thickness.
        cutIS = cutOccIS->intersectionToCut;
    }
    else
    {
        cutIS = cutOccIS->intersectionToCut->second.mirror;
    }
    double diffAngle = cutIS->second.otherTrajectory->base.region->intersection_angle(
    cutIS->second.otherTrajectory, cutIS->second.mirror->second.otherTrajectory);
    if (diffAngle > 0.5 * PI)
    {
        diffAngle = PI - diffAngle;
    }
    // don't do anything if the angle between MTs is smaller than the cutoff angle
    if (diffAngle <= (p.crossSeveringStartAngle) * PI / 180.)
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
    cutSegNumber = randomGen.rand_int(cutIS->second.occupancy - 1); // nummers 0 - occupancy -1 .
    if (cutSegNumber < cutIS->second.occupancy / 2)                // search from the beginning.
                                                    // In principle, all searches could start from the beginning. Also
                                                    // starting from the end should save (on average) half of the time
    {
        testSeg = cutIS->second.mirror->second.otherTrajectory->segments.begin();
        finalTestSeg = cutIS->second.mirror->second.otherTrajectory->segments.end();
        cutSegNumber++;
        while (testSeg != finalTestSeg) // NOTE: test should not be necessary "while(true)"
        {
            if ((**testSeg).crosses_intersection(cutIS))
            {
                cutSegNumber--;
                if (cutSegNumber == 0)
                {
                    cutSeg = *testSeg;
                    break;
                }
            }
            testSeg++;
        }
    }
    else // search from the end
    {
        testSeg = cutIS->second.mirror->second.otherTrajectory->segments.end();
        finalTestSeg = cutIS->second.mirror->second.otherTrajectory->segments.begin();
        while (testSeg != finalTestSeg)
        {
            testSeg--;
            if ((**testSeg).crosses_intersection(cutIS))
            {
                cutSegNumber++;
                if (cutSegNumber == cutIS->second.occupancy)
                {
                    cutSeg = *testSeg;
                    break;
                }
            }
        }
    }
    if (cutSeg != nullptr)
    {
        cutMT = cutSeg->mt;
        cutMT->update_length();
        if (cutMT != nullptr)
        {
            cutMT->sever_at_cross(cutIS, cutSeg);
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
    if ((p.restrictedPool != 0) && (p.vTM / (p.vPlus * (1. - totalLength / (p.poolDensity * geometry->area))) > 0.9999))
    {
        return;
    }

    // no nucleation at forbidden zone
    TrajectoryVector tv;

    // get random nucleation point (vector)
    SurfaceVector sv = geometry->random_surface_vector();

    // forbiden zone and parameter modified : no microtubule created
    if (sv.region->faceTag == -1)
    {
        return;
    }

    double overrideAngle(0.0);

    switch (p.nucleationType)
    {
            // select isotropic nucleation
        case NUC_ISOTROPIC:
            overrideAngle = sv.angle;
            break;
        default:
            break;
    }

    sv.angle = overrideAngle;
    tv = geometry->create_trajectory(sv);
    growing_mts.create(this, tv);

    return;
}

void System::determine_stochastic_event()
{
    // for(regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);

    double nextEventInterval(0.0);
    double invTotalRate(0.0);
    double typeSelector(0.0);
    double logRandom(0.0);
    double temp(0.0);

    if ((growing_mts.size() == 0) || (p.restrictedPool == 0))
    {
        double constRate(0.0);
        double slopeRate(0.0);

        constRate = p.kNuc * geometry->area + p.kRes * shrinking_mts.size() + p.kSev * totalLength;

        for (int dTag = 0; dTag <= p.faceNumber; dTag++)
        {
            constRate += p.kCat * p.RegionKcatMultiplier[dTag] * growingTipsReg[dTag];
        }

#ifdef CROSS_SEV
        constRate += p.kCross * OccupiedIntersectionList.size();
#endif

        slopeRate = p.kSev * ((p.vPlus - p.vTM) * growing_mts.size() + (p.vMin - p.vTM) * shrinking_mts.size());

        // generate a random number larger than 0 with an exponential distribution
        logRandom = -log(randomGen.rand_dbl_exc());

        if (abs(slopeRate) < ZERO_CUTOFF)
        {
            nextEventInterval = logRandom / constRate;
        }
        else
        {
            if ((temp = pow(constRate, 2) + 2 * slopeRate * logRandom) < 0)
            {

                // choose a time that is always larger than the disappearance time of all MTs
                nextEventInterval = VERY_LARGE;
            }

            else
            {
                nextEventInterval = (-constRate + sqrt(temp)) / slopeRate;
            }
        }

        // it will be negative for the case with no stochastic events
        invTotalRate = 1.0 / (constRate + slopeRate * nextEventInterval);
    }

    else
    {
        double invLMax = 1.0 / (p.poolDensity * geometry->area);
        double alpha = (p.vPlus - p.vTM) * growing_mts.size() + (p.vMin - p.vTM) * shrinking_mts.size();
        double beta = p.vPlus * growing_mts.size() * invLMax;
        double L_inf = alpha / beta;

        double baseRate = p.kNuc * geometry->area + p.kRes * shrinking_mts.size();
        for (int dTag = 0; dTag <= p.faceNumber; dTag++)
        {
            baseRate += p.kCat * p.RegionKcatMultiplier[dTag] * growingTipsReg[dTag];
        }

#ifdef CROSS_SEV
        baseRate += p.kCross * OccupiedIntersectionList.size();
#endif

        double x = p.kSev * (totalLength - L_inf) / (baseRate + p.kSev * L_inf);
        double y = -beta * log(randomGen.rand_dbl_exc()) / (baseRate + p.kSev * L_inf);

        nextEventInterval = (-x + y + lambert_w(x * exp(x - y))) / beta;

        double totalRate = baseRate + p.kSev * (L_inf + (totalLength - L_inf) * exp(-beta * nextEventInterval));
        invTotalRate = 1.0 / totalRate;
    }

#ifdef DBG_ACID_TEST
    if ((nextEventInterval < 0) || (nextEventInterval != nextEventInterval))
    {
        cerr << "Invalid deterministic event time calculated. Exiting\n";
        exit(-1);
    }
#endif

    nextStochasticEventTime = systemTime + nextEventInterval;

    // even if  invTotalRate is negative! due to the timing, the event will never be executed.
    typeSelector = randomGen();

    double typeSelectRegSource(0.0);
    for (int dTag = 0; dTag <= p.faceNumber; dTag++)
    {
        typeSelectRegSource += p.kCat * p.RegionKcatMultiplier[dTag] * growingTipsReg[dTag];
    }

    if (typeSelector < p.kNuc * geometry->area * invTotalRate)
    {
        nextStochasticEventType = NUCLEATION;
    }

    // select catastrophe event
    else if (typeSelector < (p.kNuc * geometry->area + typeSelectRegSource) * invTotalRate)
    {
        nextStochasticEventType = CATASTROPHE;
    }

    // select rescue event
    else if (typeSelector
             < (p.kNuc * geometry->area + typeSelectRegSource + p.kRes * shrinking_mts.size()) * invTotalRate)
    {
        nextStochasticEventType = RESCUE;
    }

// select severing event
#ifdef CROSS_SEV
    else if (typeSelector < (p.kNuc * geometry->area + typeSelectRegSource + p.kRes * shrinking_mts.size()
                             + p.kCross * OccupiedIntersectionList.size())
                            * invTotalRate)
    {
        nextStochasticEventType = SEVERING_AT_CROSS;
    }
#endif

    // select katanin induced event
    else
    {
        nextStochasticEventType = KATANIN;
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

        case MEASURE:
            break;
        case STATUS:
            if (!integrity_check())
            {
                cerr << "Integrity check failed. Exiting.\n";
                exit(-1);
            }

            perform_measurement();
            nextStatusEventTime = systemTime + p.measurementInterval;
            timeQueue.push_global(nextStatusEventTime, STATUS);
            break;

        case SNAPSHOT:
            // declaring variable within case: so the local blocking
            {
                OrderParameters order;

                geometry->get_order_parameters(order);

                // write a movie frame
                movieFile << order.Rdirector[0] << "," << order.Rdirector[1] << "," << order.Rdirector[2] << "#";
                movieFile << geometry->objectCM(0, 0) << "," << geometry->objectCM(1, 0) << ","
                          << geometry->objectCM(2, 0) << "#";
                movieFile << geometry->objectPA(0, 0) << "," << geometry->objectPA(1, 0) << ","
                          << geometry->objectPA(2, 0) << "*";

                geometry->output_snapshot(movieFile);

                // write a local order parameter heatmap frame
                geometry->output_order_heat_map(heatMapFile, order.localOrder, order.Sv);

                nextSnapshotEventTime = systemTime + p.movieFrameInterval;
                timeQueue.push_global(nextSnapshotEventTime, SNAPSHOT);
            }
            break;

        case STOP:
        {
            OrderParameters order;
            geometry->get_order_parameters(order);

            // division plane: polygon edge-list & area
            Eigen::Vector3d normalPPB(order.Rdirector[0], order.Rdirector[1], order.Rdirector[2]);
            geometry->areaPPB = intersecting_polygon(geometry->ppb,
                                                    geometry->globalVertex,
                                                    geometry->regions,
                                                    geometry->ppbEdgeList,
                                                    normalPPB,
                                                    geometry->nucleousPosition,
                                                    p.geometry);

            if (p.showOutput == 1)
            {
                cout << "Division plane area: " << geometry->areaPPB << endl;
            }

            stopSignal = true;
        }
        break;

        case PARAMETER_CHANGE:
        {
            cout << "Changing parameters. Reading from file [" << p.newParameterFile << "]\n";
            flush_and_reload(false);
            string oldDir = p.outputDir;
            if (p.reinitialize(p.newParameterFile.c_str()))
            {
                if (p.outputDir != oldDir)
                {
                    cout << "Switching directories.\n";
                    write_measurements_to_file(0);
                    close_files();
                    initialize_output();
                }
                p.write_to_file();

                if (p.movieEnabled)
                {
                    if (nextSnapshotEventTime < ZERO_CUTOFF)
                    {
                        timeQueue.push_global(p.movieFrameInterval, SNAPSHOT);
                        nextSnapshotEventTime = p.movieFrameInterval;
                    }
                }
                else
                {
                    nextSnapshotEventTime = 0;
                }

                if (p.newParameterReadInterval > ZERO_CUTOFF)
                {
                    timeQueue.push_global(p.newParameterReadInterval, PARAMETER_CHANGE);
                    nextParameterEventTime = p.newParameterReadInterval;
                }
            }
            else
            {
                cerr << "Failed to read parameters from file. Exiting.\n";
                emergency_break();
                exit(-1);
            }

            timeQueue.push_global(p.newParameterReadInterval, PARAMETER_CHANGE);
            nextParameterEventTime = p.newParameterReadInterval;
            break;
        }

        case SIGNAL_PPB:
            cout << " Starting formation of PPB ...\n";
            p.write_to_file();
            p.PPBkNucFraction = 0.8;
            cout << p.PPBkNucFraction << " : " << p.PPB << endl;
            OrderParameters order;
            geometry->get_order_parameters(order);
            Eigen::Vector3d normalPPB(order.Rdirector[0], order.Rdirector[1], order.Rdirector[2]);
            geometry->areaPPB = intersecting_polygon(geometry->ppb,
                                                    geometry->globalVertex,
                                                    geometry->regions,
                                                    geometry->ppbEdgeList,
                                                    normalPPB,
                                                    geometry->nucleousPosition,
                                                    p.geometry);

            for (int dTag = 0; dTag <= p.faceNumber; dTag++)
            {
                growingTipsReg[dTag] = 0;
            }

            for (size_t eno = 0; eno < geometry->regions.size(); eno++)
            {
                if (geometry->regions[eno]->polyIntersectMark == 1)
                {
                    geometry->regions[eno]->faceTag = 1;

                    for (size_t tip = 0; tip < geometry->regions[eno]->growingPlusTipList.size(); tip++)
                    {
                        growingTipsReg[1] += 1;
                    }

                    geometry->patchArea[1] += geometry->regions[eno]->area;
                    geometry->RegionsIndex[1].push_back(eno);
                }

                else
                {
                    geometry->regions[eno]->faceTag = 2;

                    for (size_t tip = 0; tip < geometry->regions[eno]->growingPlusTipList.size(); tip++)
                    {
                        growingTipsReg[2] += 1;
                    }

                    geometry->patchArea[2] += geometry->regions[eno]->area;
                    geometry->RegionsIndex[2].push_back(eno);
                }
            }
            timeQueue.push_global(VERY_LARGE, SIGNAL_PPB);
            nextPPBEventTime = VERY_LARGE;

            break;
    }

    return;
}

EventDescriptorIndex System::register_event_descriptor(EventDescriptor* ei)
{
    // maybe use a pool of numbers instead of an increasing number....this could wrap around
    EventDescriptorIndex index = EventDescriptorID++;
    EventDescriptorMap[index] = ei;

    if (EventDescriptorID == numeric_limits<EventDescriptorIndex>::max() - 10)
    {
        timeQueue.push_global(systemTime, STOP);
    }

    return index;
}

void System::unregister_event_descriptor(EventDescriptorIndex ei)
{
    EventDescriptorMap.erase(EventDescriptorMap.find(ei));
    return;
}

EventDescriptor* System::get_event_descriptor(EventDescriptorIndex idx)
{
    map<EventDescriptorIndex, EventDescriptor*>::iterator itr = EventDescriptorMap.find(idx);
    if (itr == EventDescriptorMap.end())
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
    OccupiedIntersectionList.remove_element(is.occupiedListPtr->index);
    is.occupiedListPtr = nullptr;
    is.mirror->second.occupiedListPtr = nullptr;
    return;
}

void System::add_occupied_intersection(IntersectionItr is)
{
    is->second.occupiedListPtr = OccupiedIntersectionList.create(is);
    is->second.mirror->second.occupiedListPtr = is->second.occupiedListPtr;
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
            binomialTable[i][j] = factorial[i] / (factorial[j] * factorial[i - j]);
        }
    }
    return;
}
