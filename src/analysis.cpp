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
    return o << m.m_time << "\t" << m.m_number_of_m_ts << "\t" << m.m_average_length << "\t" << m.m_length_density << "\t"
             << m.m_optical_density << "\t" << m.m_order.m_r << "\t" << m.m_order.m_rdirector[0] << "\t" << m.m_order.m_rdirector[1]
             << "\t" << m.m_order.m_rdirector[2] << "\r\n";
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
    if (m_p.m_create_subdir)
    {
        string identifier;
        time_t time_raw;
        struct tm* time_struct;
        char time_name[100];

        time(&time_raw);
        time_struct = localtime(&time_raw);
        strftime(time_name, 100, "%c", time_struct);

        // create subdirectory with a time Tag, to isolate each data folder
        m_p.m_output_dir += "/" + m_p.m_geometry + "/outputData" + string(time_name);
        int i = std::filesystem::create_directories(m_p.m_output_dir.c_str());

        // if the sub-directory path is invalid
        if (i)
        {
            cerr << "ERROR: creating output directory " << m_p.m_output_dir << " ; " << strerror(i) << "\n";
            exit(-1);
        }
    }

    // open output parameter file
    string file_name = string(m_p.m_output_dir + "/" + m_p.m_geometry + "/outputData" + string("/parameters.txt"));
    m_parameter_file.open(file_name.c_str());
    if (!m_parameter_file)
    {
        cout << "ERROR: cannot open file for output (parameters)\n";
        exit(-1);
    }

    // open output measurements file
    file_name = string(m_p.m_output_dir + "/" + m_p.m_geometry + "/outputData" + string("/measurements.txt"));
    m_measurement_file.open(file_name.c_str());
    if (!m_measurement_file.is_open())
    {
        cout << "ERROR: cannot open file for output (measurements)\n";
        exit(-1);
    }

    // write in measurement file
    write_measurement_descriptors(m_measurement_file);

    // open output Draw_MT_segments file
    file_name = string(m_p.m_output_dir + "/" + m_p.m_geometry + "/outputData" + string("/Draw_MT_segments.txt"));
    m_movie_file.open(file_name.c_str());
    if (!m_movie_file.is_open())
    {
        cout << "ERROR: cannot open file for output (Draw_MT_segments.txt)\n";
        exit(-1);
    }

    // open output heatMapOrder file
    file_name = string(m_p.m_output_dir + "/" + m_p.m_geometry + "/outputData" + string("/heatMapOrder.txt"));
    m_heat_map_file.open(file_name.c_str());
    if (!m_heat_map_file.is_open())
    {
        cout << "ERROR: cannot open file for output (heatMapOrder.txt)\n";
        exit(-1);
    }

    return;
}

void System::close_files()
{
    /// close the output files before the program terminates
    m_parameter_file.close();
    m_movie_file.close();
    m_measurement_file.close();

    return;
}

void System::perform_measurement(void)
{
    Measurement m;

    // measurementHistory reach at MAX_HISTORY_SIZE, archive only upto last MIN_HISTORY_SIZE arrays, and write down the
    // rest on output file (immediately !)
    if (m_measurement_history.size() == MAX_HISTORY_SIZE)
    {
        write_measurements_to_file(MIN_HISTORY_SIZE - 1);
    }

    m.m_time = m_system_time + m_system_time_offset;

    // length densities (per area)
    m.m_length_density = m_total_length / m_geometry->m_area;
    m.m_optical_density = m_geometry->optical_length() / m_geometry->m_area;

    // total number of microtubules
    m.m_number_of_m_ts = m_growing_mts.size() + m_shrinking_mts.size();
    m.m_growing_number = m_growing_mts.size();
    m.m_shrinking_number = m_shrinking_mts.size();

    int segment_count = 0;
    Microtubule* mt;

    mt = m_growing_mts.first();
    while (mt != nullptr)
    {
        segment_count += mt->m_segments.size();
        mt = mt->next();
    }
    mt = m_shrinking_mts.first();
    while (mt != nullptr)
    {
        segment_count += mt->m_segments.size();
        mt = mt->next();
    }

    // total number of segments
    m.m_segments = segment_count - m_boundary_crossing_count;

    // total number of trajectories (both active and inactive)
    m.m_trajectories = m_geometry->trajectory_count();

    // average length of microtubules
    if (m.m_number_of_m_ts)
    {
        m.m_average_length = m_total_length / m.m_number_of_m_ts;
        m.m_segments_per_mt = static_cast<double>(m.m_segments) / m.m_number_of_m_ts;
    }
    else
    {
        m.m_average_length = 0.;
        m.m_segments_per_mt = 1.;
    }

    // order parameter
    m_geometry->get_order_parameters(m.m_order);

    // total number zippering
    m.m_zipper_count = m_total_zipper_count;

    // total number of cross-overs
    m.m_crossover_count = m_total_crossover_count;

    // total number of induced catastrophes
    m.m_induced_catastrophe_count = m_total_induced_catastrophe_count;

    // total number of valid events
    m.m_valid_d_event_count = m_total_valid_d_event_count;

    // total number of invalid events
    m.m_invalid_d_event_count = m_total_invalid_d_event_count;

    // total number of severing event
    m.m_s_event_count = m_total_s_event_count;
    m.m_length_severing_count = m_total_length_severing_count;
    m.m_intersection_severing_count = m_total_intersection_severing_count;

// total number of occupied intersections
#ifdef CROSS_SEV
    m.m_occupied_intersection_count = m_occupied_intersection_list.size();
#else
    m.occupiedIntersectionCount = 0;
#endif

    // G-effective at difference faces
    // for(int dTag =0;dTag <=p.faceNumber; dTag++)
    // m.G_effAdjust.push_back = (p.kRes/(-p.vMin+p.vTM) - p.RegionKcatMultiplier[dTag]*p.kCat/(realVplus -
    // p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/
    // (p.kNuc*(-p.vMin+realVplus)*realVplus), 1./3.);

    // archive (store) measurement
    m_measurement_history.push_back(m);

    // show output
    if (m_p.m_show_output == 1)
    {
        cout << "[t=" << setprecision(9) << m.m_time << setprecision(6) << "\t\t" << m.m_trajectories << " trajectories \t"
             << m.m_segments << " segments]\n"
             << "density [real]: " << m.m_length_density << "\t\t [optical]: " << m.m_optical_density << "\n"
             << "length/MT     : " << m.m_average_length << "\tMT#: " << m.m_number_of_m_ts
             << "\tsegments/MT: " << m.m_segments_per_mt << "\n"
             << "R  [real]     : " << m.m_order.m_r << "[" << m.m_order.m_rdirector[0] << ", " << m.m_order.m_rdirector[1] << ", "
             << m.m_order.m_rdirector[2] << "]" << "\n"
             << "C  [real]     : " << m.m_order.m_c << "\n\n";
        //<< "G_eff_adjust  : " << m.G_effAdjust[0] << " [on regular surfaces]\n\n";
    }

    return;
}

void System::write_measurements_to_file(int number_to_keep)
{
    int i, i_max;

    i_max = m_measurement_history.size() - number_to_keep;

    // write output on files
    for (i = 0; i < i_max; i++)
    {
        m_measurement_file << m_measurement_history[0];
        m_measurement_history.pop_front();
    }

    return;
}
