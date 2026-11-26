#include <filesystem>
#include "types.hpp"
#include "measurement.hpp"
#include "geometry.hpp"
#include "system.hpp"
#include "linalg.hpp"

ostream& operator<<(ostream& o, const Measurement m)
{
    // set precession to 4 decimal digits
    o.setf(ios::fixed, ios::floatfield);
    o.precision(4);

    // write output data in file
    return o << m.time << "\t" << m.number_of_m_ts << "\t" << m.average_length << "\t" << m.length_density << "\t"
             << m.optical_density << "\t" << m.order.r << "\t" << m.order.rdirector[0] << "\t" << m.order.rdirector[1]
             << "\t" << m.order.rdirector[2] << "\r\n";
}

void write_measurement_descriptors(ostream& o)
{
    // write header of the output file
    o << "time\t\t" << " #MT\t" << " <l>\t" << " rho\t" << "rho_opt\t" << " R\t" << "R_x\t" << "R_y\t" << "R_z\r\n";
    return;
}

void System::initialize_output(void)
{

    // if create sub-directory in enabled
    if (p.create_subdir)
    {
        string identifier;
        time_t time_raw;
        struct tm* time_struct;
        char time_name[100];

        time(&time_raw);
        time_struct = localtime(&time_raw);
        strftime(time_name, 100, "%c", time_struct);

        // create subdirectory with a time Tag, to isolate each data folder
        p.output_dir += "/" + p.geometry + "/outputData" + string(time_name);
        int i = std::filesystem::create_directories(p.output_dir.c_str());

        // if the sub-directory path is invalid
        if (i)
        {
            cerr << "ERROR: creating output directory " << p.output_dir << " ; " << strerror(i) << "\n";
            exit(-1);
        }
    }

    // open output parameter file
    string file_name = string(p.output_dir + "/" + p.geometry + "/outputData" + string("/parameters.txt"));
    parameter_file.open(file_name.c_str());
    if (!parameter_file)
    {
        cout << "ERROR: cannot open file for output (parameters)\n";
        exit(-1);
    }

    // open output measurements file
    file_name = string(p.output_dir + "/" + p.geometry + "/outputData" + string("/measurements.txt"));
    measurement_file.open(file_name.c_str());
    if (!measurement_file.is_open())
    {
        cout << "ERROR: cannot open file for output (measurements)\n";
        exit(-1);
    }

    // write in measurement file
    write_measurement_descriptors(measurement_file);

    // open output Draw_MT_segments file
    file_name = string(p.output_dir + "/" + p.geometry + "/outputData" + string("/Draw_MT_segments.txt"));
    movie_file.open(file_name.c_str());
    if (!movie_file.is_open())
    {
        cout << "ERROR: cannot open file for output (Draw_MT_segments.txt)\n";
        exit(-1);
    }

    // open output heatMapOrder file
    file_name = string(p.output_dir + "/" + p.geometry + "/outputData" + string("/heatMapOrder.txt"));
    heat_map_file.open(file_name.c_str());
    if (!heat_map_file.is_open())
    {
        cout << "ERROR: cannot open file for output (heatMapOrder.txt)\n";
        exit(-1);
    }

    return;
}

void System::close_files()
{
    /// close the output files before the program terminates
    parameter_file.close();
    movie_file.close();
    measurement_file.close();

    return;
}

void System::perform_measurement(void)
{
    Measurement m;

    // measurementHistory reach at MAX_HISTORY_SIZE, archive only upto last MIN_HISTORY_SIZE arrays, and write down the
    // rest on output file (immediately !)
    if (measurement_history.size() == MAX_HISTORY_SIZE)
    {
        write_measurements_to_file(MIN_HISTORY_SIZE - 1);
    }

    m.time = system_time + system_time_offset;

    // length densities (per area)
    m.length_density = total_length / geometry->area;
    m.optical_density = geometry->optical_length() / geometry->area;

    // total number of microtubules
    m.number_of_m_ts = growing_mts.size() + shrinking_mts.size();
    m.growing_number = growing_mts.size();
    m.shrinking_number = shrinking_mts.size();

    int segment_count = 0;
    Microtubule* mt;

    mt = growing_mts.first();
    while (mt != nullptr)
    {
        segment_count += mt->segments.size();
        mt = mt->next();
    }
    mt = shrinking_mts.first();
    while (mt != nullptr)
    {
        segment_count += mt->segments.size();
        mt = mt->next();
    }

    // total number of segments
    m.segments = segment_count - boundary_crossing_count;

    // total number of trajectories (both active and inactive)
    m.trajectories = geometry->trajectory_count();

    // average length of microtubules
    if (m.number_of_m_ts)
    {
        m.average_length = total_length / m.number_of_m_ts;
        m.segments_per_mt = static_cast<double>(m.segments) / m.number_of_m_ts;
    }
    else
    {
        m.average_length = 0.;
        m.segments_per_mt = 1.;
    }

    // order parameter
    geometry->get_order_parameters(m.order);

    // total number zippering
    m.zipper_count = total_zipper_count;

    // total number of cross-overs
    m.crossover_count = total_crossover_count;

    // total number of induced catastrophes
    m.induced_catastrophe_count = total_induced_catastrophe_count;

    // total number of valid events
    m.valid_d_event_count = total_valid_d_event_count;

    // total number of invalid events
    m.invalid_d_event_count = total_invalid_d_event_count;

    // total number of severing event
    m.s_event_count = total_s_event_count;
    m.length_severing_count = total_length_severing_count;
    m.intersection_severing_count = total_intersection_severing_count;

// total number of occupied intersections
#ifdef CROSS_SEV
    m.occupied_intersection_count = occupied_intersection_list.size();
#else
    m.occupiedIntersectionCount = 0;
#endif

    // G-effective at difference faces
    // for(int dTag =0;dTag <=p.faceNumber; dTag++)
    // m.G_effAdjust.push_back = (p.kRes/(-p.vMin+p.vTM) - p.RegionKcatMultiplier[dTag]*p.kCat/(realVplus -
    // p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/
    // (p.kNuc*(-p.vMin+realVplus)*realVplus), 1./3.);

    // archive (store) measurement
    measurement_history.push_back(m);

    // show output
    if (p.show_output == 1)
    {
        cout << "[t=" << setprecision(9) << m.time << setprecision(6) << "\t\t" << m.trajectories << " trajectories \t"
             << m.segments << " segments]\n"
             << "density [real]: " << m.length_density << "\t\t [optical]: " << m.optical_density << "\n"
             << "length/MT     : " << m.average_length << "\tMT#: " << m.number_of_m_ts
             << "\tsegments/MT: " << m.segments_per_mt << "\n"
             << "R  [real]     : " << m.order.r << "[" << m.order.rdirector[0] << ", " << m.order.rdirector[1] << ", "
             << m.order.rdirector[2] << "]" << "\n"
             << "C  [real]     : " << m.order.c << "\n\n";
        //<< "G_eff_adjust  : " << m.G_effAdjust[0] << " [on regular surfaces]\n\n";
    }

    return;
}

void System::write_measurements_to_file(int number_to_keep)
{
    int i, i_max;

    i_max = measurement_history.size() - number_to_keep;

    // write output on files
    for (i = 0; i < i_max; i++)
    {
        measurement_file << measurement_history[0];
        measurement_history.pop_front();
    }

    return;
}
