#include "corticalSimReal.h"

ostream& operator<<(ostream& o, const Measurement m)
{
    // set precession to 4 decimal digits
    o.setf(ios::fixed,ios::floatfield);
    o.precision(4);

    // write output data in file
    return o << m.time << "\t"<< m.numberOfMTs <<"\t"<< m.averageLength << "\t"<< m.lengthDensity << "\t"
             << m.opticalDensity << "\t"<< m.order.R << "\t" << m.order.Rdirector[0] << "\t"
             << m.order.Rdirector[1] << "\t" << m.order.Rdirector[2] <<"\r\n";

}

void writeMeasurementDescriptors(ostream& o)
{
    // write header of the output file
    o << "time\t\t" <<" #MT\t"<< " <l>\t"<<" rho\t"<<"rho_opt\t"<<" R\t"<<"R_x\t"<<"R_y\t"<<"R_z\r\n";
    return;
}

void System::initializeOutput(void)
{
    int j;

    // if create sub-directory in enabled
    if (p.createSubdir)
    {
        string identifier;
        time_t timeRaw;
        struct tm * timeStruct;
        char timeName[100];

        time(&timeRaw);
        timeStruct = localtime(&timeRaw);
        strftime(timeName,100,"%c",timeStruct);

        // create subdirectory with a time Tag, to isolate each data folder
        p.outputDir += "/" + p.geometry + "/outputData" + string(timeName);
        int i = boost::filesystem::create_directories(p.outputDir.c_str());

        // if the sub-directory path is invalid
        if (i)
        {
            cerr << "ERROR: creating output directory " << p.outputDir << " ; " << strerror(i) << "\n";
            exit(-1);
        }
    }

    // open output parameter file
    string FileName = string(p.outputDir + "/" + p.geometry + "/outputData" + string("/parameters.txt"));
    parameterFile.open(FileName.c_str());
    if (!parameterFile)
    {
        cout << "ERROR: cannot open file for output (parameters)\n";
        exit(-1);
    }

    // open output measurements file
    FileName = string(p.outputDir + "/" + p.geometry + "/outputData" +  string("/measurements.txt"));
    measurementFile.open(FileName.c_str());
    if (!measurementFile.is_open())
    {
        cout << "ERROR: cannot open file for output (measurements)\n";
        exit(-1);
    }

    // write in measurement file
    writeMeasurementDescriptors(measurementFile);

    // open output Draw_MT_segments file
    FileName = string(p.outputDir + "/" + p.geometry + "/outputData" +  string("/Draw_MT_segments.txt"));
    movieFile.open(FileName.c_str());
    if (!movieFile.is_open())
    {
        cout << "ERROR: cannot open file for output (Draw_MT_segments.txt)\n";
        exit(-1);
    }

    // open output heatMapOrder file
    FileName = string(p.outputDir + "/" + p.geometry + "/outputData" + string("/heatMapOrder.txt"));
    heatMapFile.open(FileName.c_str());
    if (!heatMapFile.is_open())
    {
        cout << "ERROR: cannot open file for output (heatMapOrder.txt)\n";
        exit(-1);
    }

    return;
}

void System::closeFiles()
{
    /// close the output files before the program terminates
    parameterFile.close();
    movieFile.close();
    measurementFile.close();

    return;
}

void System::performMeasurement(void)
{
    Measurement m;
    double *h;
    double *al;
    int *an;
    int i;

    // measurementHistory reach at MAX_HISTORY_SIZE, archive only upto last MIN_HISTORY_SIZE arrays, and write down the rest on output file (immediately !)
    if (measurementHistory.size() == MAX_HISTORY_SIZE)
    {
        writeMeasurementsToFile(MIN_HISTORY_SIZE-1);
    }

    m.time = systemTime + systemTimeOffset;

    // length densities (per area)
    m.lengthDensity = totalLength/geometry->area;
    m.opticalDensity = geometry->opticalLength()/geometry->area;

    // total number of microtubules
    m.numberOfMTs = growing_mts.size() + shrinking_mts.size();
    m.growingNumber = growing_mts.size();
    m.shrinkingNumber = shrinking_mts.size();

    int segmentCount = 0;
    Microtubule* mt;

    mt = growing_mts.first();
    while (mt != NULL)
    {
        segmentCount += mt->segments.size();
        mt = mt->next();
    }
    mt = shrinking_mts.first();
    while (mt != NULL)
    {
        segmentCount += mt->segments.size();
        mt = mt->next();
    }

    // total number of segments
    m.segments = segmentCount - boundaryCrossingCount;

    // total number of trajectories (both active and inactive)
    m.trajectories = geometry->trajectoryCount();

    // average length of microtubules
    if(m.numberOfMTs)
    {
        m.averageLength = totalLength/m.numberOfMTs;
        m.segmentsPerMT = static_cast<double>(m.segments)/m.numberOfMTs;
    }
    else
    {
        m.averageLength = 0.;
        m.segmentsPerMT = 1.;
    }

    // order parameter
    geometry->getOrderParameters(m.order);

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

    // value of plus end velocity
    double realVplus;

    // infinite tubulin pool
    if (p.restrictedPool == 0)
        realVplus = p.vPlus;
    // finite tubulin pool
    else
        realVplus = p.vPlus*(1.0 - totalLength/(p.poolDensity*geometry->area));

    // G-effective at difference faces
    //for(int dTag =0;dTag <=p.faceNumber; dTag++)
    //m.G_effAdjust.push_back = (p.kRes/(-p.vMin+p.vTM) - p.RegionKcatMultiplier[dTag]*p.kCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/
			  // (p.kNuc*(-p.vMin+realVplus)*realVplus), 1./3.);

    // archive (store) measurement
    measurementHistory.push_back(m);

   // show output
   if(p.showOutput ==1)
   {
         cout << "[t=" << setprecision(9) << m.time << setprecision(6) << "\t\t" << m.trajectories << " trajectories \t" << m.segments << " segments]\n"
         << "density [real]: " << m.lengthDensity << "\t\t [optical]: " << m.opticalDensity << "\n"
         << "length/MT     : " << m.averageLength << "\tMT#: " << m.numberOfMTs << "\tsegments/MT: " << m.segmentsPerMT << "\n"
         << "R  [real]     : " << m.order.R <<"[" << m.order.Rdirector[0] << ", " << m.order.Rdirector[1] << ", " << m.order.Rdirector[2] << "]"<<"\n"
         << "C  [real]     : " << m.order.C <<"\n\n";
         //<< "G_eff_adjust  : " << m.G_effAdjust[0] << " [on regular surfaces]\n\n";
    }

    return;
}

void System::writeMeasurementsToFile(int numberToKeep)
{
    int i,j, iMax;

    iMax = measurementHistory.size()-numberToKeep;

    // write output on files
    for (i=0; i<iMax; i++)
    {
        measurementFile << measurementHistory[0];
        measurementHistory.pop_front();
    }

    return;
}

