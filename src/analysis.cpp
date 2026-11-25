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
    return o << m.time << "\t" << m.numberOfMTs << "\t" << m.averageLength << "\t" << m.lengthDensity << "\t"
             << m.opticalDensity << "\t" << m.order.R << "\t" << m.order.Rdirector[0] << "\t" << m.order.Rdirector[1]
             << "\t" << m.order.Rdirector[2] << "\r\n";
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
    if (p.createSubdir)
    {
        string identifier;
        time_t time_raw;
        struct tm* time_struct;
        char time_name[100];

        time(&time_raw);
        time_struct = localtime(&time_raw);
        strftime(time_name, 100, "%c", time_struct);

        // create subdirectory with a time Tag, to isolate each data folder
        p.outputDir += "/" + p.geometry + "/outputData" + string(time_name);
        int i = std::filesystem::create_directories(p.outputDir.c_str());

        // if the sub-directory path is invalid
        if (i)
        {
            cerr << "ERROR: creating output directory " << p.outputDir << " ; " << strerror(i) << "\n";
            exit(-1);
        }
    }

    // open output parameter file
    string file_name = string(p.outputDir + "/" + p.geometry + "/outputData" + string("/parameters.txt"));
    parameterFile.open(file_name.c_str());
    if (!parameterFile)
    {
        cout << "ERROR: cannot open file for output (parameters)\n";
        exit(-1);
    }

    // open output measurements file
    file_name = string(p.outputDir + "/" + p.geometry + "/outputData" + string("/measurements.txt"));
    measurementFile.open(file_name.c_str());
    if (!measurementFile.is_open())
    {
        cout << "ERROR: cannot open file for output (measurements)\n";
        exit(-1);
    }

    // write in measurement file
    write_measurement_descriptors(measurementFile);

    // open output Draw_MT_segments file
    file_name = string(p.outputDir + "/" + p.geometry + "/outputData" + string("/Draw_MT_segments.txt"));
    movieFile.open(file_name.c_str());
    if (!movieFile.is_open())
    {
        cout << "ERROR: cannot open file for output (Draw_MT_segments.txt)\n";
        exit(-1);
    }

    // open output heatMapOrder file
    file_name = string(p.outputDir + "/" + p.geometry + "/outputData" + string("/heatMapOrder.txt"));
    heatMapFile.open(file_name.c_str());
    if (!heatMapFile.is_open())
    {
        cout << "ERROR: cannot open file for output (heatMapOrder.txt)\n";
        exit(-1);
    }

    return;
}

void System::close_files()
{
    /// close the output files before the program terminates
    parameterFile.close();
    movieFile.close();
    measurementFile.close();

    return;
}

void System::perform_measurement(void)
{
    Measurement m;

    // measurementHistory reach at MAX_HISTORY_SIZE, archive only upto last MIN_HISTORY_SIZE arrays, and write down the
    // rest on output file (immediately !)
    if (measurementHistory.size() == MAX_HISTORY_SIZE)
    {
        write_measurements_to_file(MIN_HISTORY_SIZE - 1);
    }

    m.time = systemTime + systemTimeOffset;

    // length densities (per area)
    m.lengthDensity = totalLength / geometry->area;
    m.opticalDensity = geometry->optical_length() / geometry->area;

    // total number of microtubules
    m.numberOfMTs = growing_mts.size() + shrinking_mts.size();
    m.growingNumber = growing_mts.size();
    m.shrinkingNumber = shrinking_mts.size();

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
    m.segments = segment_count - boundaryCrossingCount;

    // total number of trajectories (both active and inactive)
    m.trajectories = geometry->trajectory_count();

    // average length of microtubules
    if (m.numberOfMTs)
    {
        m.averageLength = totalLength / m.numberOfMTs;
        m.segmentsPerMT = static_cast<double>(m.segments) / m.numberOfMTs;
    }
    else
    {
        m.averageLength = 0.;
        m.segmentsPerMT = 1.;
    }

    // order parameter
    geometry->get_order_parameters(m.order);

    // total number zippering
    m.zipperCount = totalZipperCount;

    // total number of cross-overs
    m.crossoverCount = totalCrossoverCount;

    // total number of induced catastrophes
    m.inducedCatastropheCount = totalInducedCatastropheCount;

    // total number of valid events
    m.validDEventCount = totalValidDEventCount;

    // total number of invalid events
    m.invalidDEventCount = totalInvalidDEventCount;

    // total number of severing event
    m.sEventCount = totalSEventCount;
    m.lengthSeveringCount = totalLengthSeveringCount;
    m.intersectionSeveringCount = totalIntersectionSeveringCount;

// total number of occupied intersections
#ifdef CROSS_SEV
    m.occupiedIntersectionCount = OccupiedIntersectionList.size();
#else
    m.occupiedIntersectionCount = 0;
#endif

    // G-effective at difference faces
    // for(int dTag =0;dTag <=p.faceNumber; dTag++)
    // m.G_effAdjust.push_back = (p.kRes/(-p.vMin+p.vTM) - p.RegionKcatMultiplier[dTag]*p.kCat/(realVplus -
    // p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/
    // (p.kNuc*(-p.vMin+realVplus)*realVplus), 1./3.);

    // archive (store) measurement
    measurementHistory.push_back(m);

    // show output
    if (p.showOutput == 1)
    {
        cout << "[t=" << setprecision(9) << m.time << setprecision(6) << "\t\t" << m.trajectories << " trajectories \t"
             << m.segments << " segments]\n"
             << "density [real]: " << m.lengthDensity << "\t\t [optical]: " << m.opticalDensity << "\n"
             << "length/MT     : " << m.averageLength << "\tMT#: " << m.numberOfMTs
             << "\tsegments/MT: " << m.segmentsPerMT << "\n"
             << "R  [real]     : " << m.order.R << "[" << m.order.Rdirector[0] << ", " << m.order.Rdirector[1] << ", "
             << m.order.Rdirector[2] << "]" << "\n"
             << "C  [real]     : " << m.order.C << "\n\n";
        //<< "G_eff_adjust  : " << m.G_effAdjust[0] << " [on regular surfaces]\n\n";
    }

    return;
}

void System::write_measurements_to_file(int number_to_keep)
{
    int i, i_max;

    i_max = measurementHistory.size() - number_to_keep;

    // write output on files
    for (i = 0; i < i_max; i++)
    {
        measurementFile << measurementHistory[0];
        measurementHistory.pop_front();
    }

    return;
}
