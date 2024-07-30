#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <iomanip>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "corticalSim.h"


ostream& operator<<(ostream& o, const Measurement m)
{
#ifdef BAND_CAT
	return o << setprecision(9) << m.time << setprecision(6) << "\t" << m.lengthDensity << "\t" 	\
				<< m.averageLength << "\t" 	\
				<< m.order.S2 << "\t" << m.order.S2angle << "\t" << m.order.S4 << "\t" << m.order.S4angle << "\t" \
				<< m.growingNumber << "\t" << m.shrinkingNumber << "\t" << m.segments << "\t" << m.trajectories << "\t" \
				<< m.zipperCount << "\t" << m.crossoverCount << "\t" << m.inducedCatastropheCount << "\t"\
				<< m.validDEventCount << "\t" << m.invalidDEventCount << "\t" << m.sEventCount  << "\t" \
				<< m.opticalDensity << "\t"  << m.numberOfMTs << "\t" << m.segmentsPerMT<<  "\t" << m.lengthSeveringCount << "\t" \
				<< m.intersectionSeveringCount << "\t" << m.order.S2Opt << "\t" << m.order.S2angleOpt << "\t" \
				<< m.order.S4Opt << "\t" << m.order.S4angleOpt << "\t" << m.occupiedIntersectionCount << "\t" \
				<< m.order.R << "\t" << m.order.Rdirector[0] << "\t" << m.order.Rdirector[1] << "\t" << m.order.Rdirector[2] << "\t" \
				<< m.G_effAdjust_normal << "\t" << m.G_effAdjust_special << "\t" << m.G_effMeasured << "\t" << m.G_effAdjust_band << "\t" << m.G_effAdjust_gap << "\t" \
				<< m.G_effAdjust_min << "\t" << m.G_effAdjust_max << "\t"	<< m.G_effAdjust_S2angle << "\t"	<< m.nucleationCount  << "\t"	<< m.deflectionCount << "\t"	<< m.nOccupiedNCs << "\r\n";
				//<< m.G_effAdjust_normal << "\t" << m.G_effAdjust_special << "\t"	<< m.G_effAdjust_band << "\t" << m.G_effAdjust_gap << "\r\n";
#else
	return o << setprecision(9) << m.time << setprecision(6) << "\t" << m.lengthDensity << "\t" 	\
				<< m.averageLength << "\t" 	\
				<< m.order.S2 << "\t" << m.order.S2angle << "\t" << m.order.S4 << "\t" << m.order.S4angle << "\t" \
				<< m.growingNumber << "\t" << m.shrinkingNumber << "\t" << m.segments << "\t" << m.trajectories << "\t" \
				<< m.zipperCount << "\t" << m.crossoverCount << "\t" << m.inducedCatastropheCount << "\t"\
				<< m.validDEventCount << "\t" << m.invalidDEventCount << "\t" << m.sEventCount  << "\t" \
				<< m.opticalDensity << "\t"  << m.numberOfMTs << "\t" << m.segmentsPerMT<<  "\t" << m.lengthSeveringCount << "\t" \
				<< m.intersectionSeveringCount << "\t" << m.order.S2Opt << "\t" << m.order.S2angleOpt << "\t" \
				<< m.order.S4Opt << "\t" << m.order.S4angleOpt << "\t" << m.occupiedIntersectionCount << "\t" \
				<< m.order.R << "\t" << m.order.Rdirector[0] << "\t" << m.order.Rdirector[1] << "\t" << m.order.Rdirector[2] << "\t" \
				<< m.G_effAdjust_normal << "\t" << m.G_effAdjust_special << "\t" << m.G_effMeasured << "\t" \
				<< m.G_effAdjust_min << "\t" << m.G_effAdjust_max << "\t"	<< m.G_effAdjust_S2angle << "\t"	<< m.nucleationCount  << "\t"	<< m.deflectionCount << "\t"	<< m.nOccupiedNCs << "\r\n";
				//<< m.G_effAdjust_normal << "\t" << m.G_effAdjust_special << "\r\n";
#endif
}

void writeMeasurementDescriptors(ostream& o)
{
#ifdef BAND_CAT
	o << "time\tdensity"\
		<< "\t<l>" \
		<< "\tS2\tS2angle\tS4\tS4angle\t#growing\t#shrinking\t#segments\t#trajectories"\
		<< "\tzippering_ev\tcrossover_ev\tinduced_catastrophe_ev\t"\
		<< "valid_deterministic_ev\tinvalid_deterministic_ev\tstochastic_ev"\
		<< "\toptical_density\tmicrotubules\tsegments_per_MT\trandom_severing_ev\tintersection_severing_ev"\
		<< "\tS2opt\tS2opt_angle\tS4opt\tS4opt_angle\toccupied_intersections"\
		<< "\tR\tR_x\tR_y\tR_z"\
		<< "\tG_eff_adjusted_normal\tG_eff_adjusted_special\tG_eff_measured\tG_eff_adjusted_band\tG_eff_adjusted_gap\t"\
		<< "\tG_eff_adjusted_min\tG_eff_adjusted_max\tG_eff_adjusted_S2angle\tnucleation_ev\tdeflection_ev\tn_occupied_NCs\r\n";
		//<< "\tG_eff_adjusted normal\tG_eff_adjusted special\tG_eff_adjusted band\tG_eff_adjusted gap\r\n";
#else
	o << "time\tdensity"\
		<< "\t<l>" \
		<< "\tS2\tS2angle\tS4\tS4angle\t#growing\t#shrinking\t#segments\t#trajectories"\
		<< "\tzippering_ev\tcrossover_ev\tinduced_catastrophe_ev\t"\
		<< "valid_deterministic_ev\tinvalid_deterministic_ev\tstochastic_ev"\
		<< "\toptical_density\tmicrotubules\tsegments_per_MT\trandom_severing_ev\tintersection_severing_ev"\
		<< "\tS2opt\tS2opt_angle\tS4opt\tS4opt_angle\toccupied_intersections"\
		<< "\tR\tR_x\tR_y\tR_z"\
		<< "\tG_eff_adjusted_normal\tG_eff_adjusted_special\tG_eff_measured\t"\
		<< "\tG_eff_adjusted_min\tG_eff_adjusted_max\tG_eff_adjusted_S2angle\tnucleation_ev\tdeflection_ev\tn_occupied_NCs\r\n";
		//<< "\tG_eff_adjusted normal\tG_eff_adjusted special\r\n";
#endif
	return;
}

void writeNucleationPositions(string outputdir, vector<double> &nucleationpositions, double time)
{
  string filename = "/nucleationsXpos.txt";
  ofstream outfile;
  outfile.open(outputdir + filename,ios::out | ios::app);
  outfile << "time " << time << '\n';
  for (double elem : nucleationpositions)
  {
    outfile << elem << '\t';
  }
  outfile << "\n\n";
  outfile.close();
  return;
}

void writeRegionOrder(vector<OrderParameters> &ops, string outputdir, const int number, double time)
{
    vector<string> filenames = {"/localS2.txt","/localS2angle.txt","/localS2Opt.txt","/localS2angleOpt.txt","/localL.txt","/localLOpt.txt"};
    for (int fi = 0; fi < filenames.size(); fi++)
    {
        ofstream outfile;
        outfile.open(outputdir + filenames.at(fi),ios::out | ios::app);
        outfile << "time " << time << '\n';
        int row = 0;
        for(int i=0; i < ops.size(); i++)
        {
            int regionColumn = (i)%number;
            int regionRow = (i)/number;
            if (regionRow > row)
            {
                row++;
                outfile << '\n';
            }
            switch (fi)
            {
                case 0:
                    outfile << ops.at(i).S2 << '\t';
                    break;
                case 1:
                    outfile << ops.at(i).S2angle << '\t';
                    break;
                case 2:
                    outfile << ops.at(i).S2Opt << '\t';
                    break;
                case 3:
                    outfile << ops.at(i).S2angleOpt << '\t';
                    break;
                case 4:
                    outfile << ops.at(i).localL << '\t';
                    break;
                case 5:
                    outfile << ops.at(i).localLOpt << '\t';
                    break;
            }
        }
        outfile << "\n\n";
        outfile.close();
    }
    return;
}



void System::initializeOutput(void)
{
	int j;

	if (p.createSubdir)
	{
		string identifier;
		time_t timeRaw;
		struct tm * timeStruct;
		char timeName[100];	// way too much, but ok
		
		time(&timeRaw);
		timeStruct = localtime(&timeRaw);
		strftime(timeName,100,"%y%m%d-%H%M%S",timeStruct);
	
		#ifdef _WIN32
		p.outputDir += "\\" + string(timeName);
		int i = mkdir(p.outputDir.c_str());
		#else
		p.outputDir += "/" + string(timeName);
		int i = mkdir(p.outputDir.c_str(), 0777);
		#endif
		if (i)
		{
			cerr << "ERROR: creating output directory " << p.outputDir << " ; " << strerror(i) << "\n";
			exit(-1);
		}
	}

	parameterFile.open(string(p.outputDir + string("/parameters.txt")).c_str());
	if (!parameterFile)
	{
		cout << "ERROR: cannot open file for output [" << string(p.outputDir +  string("/parameters.txt")) << "]\n";
		exit(-1);
	}

	mtLengthHistogram.setFileName(p.outputDir + string("/mtLengthHistogram.txt"));
	mtLifetimeHistogram.setFileName(p.outputDir + string("/mtLifetimeHistogram.txt"));
	segAngleLengthHistogram.setFileName(p.outputDir + string("/segAngleLengthHistogram.txt"));
	segAngleLifetimeHistogram.setFileName(p.outputDir + string("/segAngleLifetimeHistogram.txt"));
    if (p.spatialHistogramType == sh_x || p.spatialHistogramType == sh_xy)
    {
        histX->setFileName(p.outputDir + string("/spatialHistogramX.txt"));
        histX->writeHeader();
    }
    if (p.spatialHistogramType == sh_y || p.spatialHistogramType == sh_xy)
    {
        histY->setFileName(p.outputDir + string("/spatialHistogramY.txt"));
        histY->writeHeader();
    }

	measurementFile.open(string(p.outputDir + string("/measurements.txt")).c_str());
	if (!measurementFile.is_open())
	{
		cout << "ERROR: cannot open file for output (measurements)\n";
		exit(-1);
	}
	writeMeasurementDescriptors(measurementFile);

	movieFile.open(string(p.outputDir + "/" + string("snapshots.txt")).c_str());
	if (!movieFile.is_open())
	{
		cout << "ERROR: cannot open file for output (snapshots)\n";
		exit(-1);
	}

	angleHistogramFile.open(string(p.outputDir + string("/angleHistogramHalf.txt")).c_str());
    if (!angleHistogramFile.is_open())
  	{
		cout << "ERROR: cannot open file for output (angle histogram)\n";
        exit(-1);
	}
    
    angleHistogramOpticalFile.open(string(p.outputDir + string("/angleHistogramHalfOptical.txt")).c_str());
    if (!angleHistogramOpticalFile.is_open())
    {
      cout << "ERROR: cannot open file for output (angle histogram optical)\n";
      exit(-1);
    }

    if ((p.discreteAngleNumber) && ((geometry->type == g_periodic) || (geometry->type == g_grid)))
    {
	 	angleLengthFile.open(string(p.outputDir + string("/discreteAngleLengths.txt")).c_str());
		angleNumberFile.open(string(p.outputDir + string("/discreteAngleNumbers.txt")).c_str());
	 
        if (!angleLengthFile.is_open() || !angleNumberFile.is_open())
	    {
	        cout << "ERROR: cannot open file for output (discrete angle lengths)\n";
	        exit(-1);
	    }
		angleLengthFile << "# Average length at angles:\t";
		for (j=0; j<p.discreteAngleNumber; j++)
        {
           angleLengthFile << *(p.nucleationAngles+j) << "\t";
        }
        angleLengthFile << "\n";
		angleNumberFile << "# Average number at angles:\t";
		for (j=0; j<p.discreteAngleNumber; j++)
        {
            angleNumberFile << *(p.nucleationAngles+j) << "\t";
        }
        angleNumberFile << "\n";

    }
    
    #ifdef MTBASED_NUCLEATION_PROBABILITY
    	MTbasedProbDensityFile.open(string(p.outputDir + string("/densityMTbased.txt")).c_str());
    	if (!MTbasedProbDensityFile.is_open())
    	{
      		cout << "ERROR: cannot open file for output (MT-based probability vs density)\n";
      		exit(-1);
    	}
    	unboundProbDensityFile.open(string(p.outputDir + string("/densityUnbound.txt")).c_str());
    	if (!unboundProbDensityFile.is_open())
    	{
      		cout << "ERROR: cannot open file for output (unbound probability vs density)\n";
      		exit(-1);
    	}
    #endif
	
	return;	
}


void System::closeFiles()
{
	// to make sure that buffers get flushed before the program is killed, close the files
	parameterFile.close();

	movieFile.close();
	measurementFile.close();

	angleHistogramFile.close();
	angleHistogramOpticalFile.close();
	angleLengthFile.close();
	angleNumberFile.close();
	#ifdef MTBASED_NUCLEATION_PROBABILITY
	MTbasedProbDensityFile.close();
	unboundProbDensityFile.close();
	#endif

	mtLengthHistogram.closeFile();
	mtLifetimeHistogram.closeFile();
	segAngleLengthHistogram.closeFile();
	segAngleLifetimeHistogram.closeFile();
    if (p.spatialHistogramType == sh_x || p.spatialHistogramType == sh_xy)
        histX->closeFile();
    if (p.spatialHistogramType == sh_y || p.spatialHistogramType == sh_xy)
        histY->closeFile();

	return;
}

void System::performMeasurement(void)
{
	Measurement m;
	double *h;
	double *al;
	int *an;
	int i;
	
	if (measurementHistory.size() == MAX_HISTORY_SIZE)
	{
		writeMeasurementsToFile(MIN_HISTORY_SIZE-1);
	}

	m.time = systemTime + systemTimeOffset;

	m.lengthDensity = totalLength/geometry->area;
	m.opticalDensity = geometry->opticalLength()/geometry->area;

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
	m.segments = segmentCount - boundaryCrossingCount;
	m.trajectories = geometry->trajectoryCount(); // NOTE: both active and inactive
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

	geometry->getOrderParameters(m.order);
  if (p.regionalOutputQuantities)
    geometry->getLocalOrderParameters(m.local_order);// for gridcylinder only


	m.zipperCount = totalZipperCount;
	m.crossoverCount = totalCrossoverCount;
	m.inducedCatastropheCount = totalInducedCatastropheCount;
	m.validDEventCount = totalValidDEventCount;
	m.invalidDEventCount = totalInvalidDEventCount;
	m.sEventCount = totalSEventCount;
	m.lengthSeveringCount = totalLengthSeveringCount;
	m.intersectionSeveringCount = totalIntersectionSeveringCount;
	m.occupiedIntersectionCount = OccupiedIntersectionList.size();
  m.nucleationCount = totalNucleationCount;
  m.deflectionCount = totalDeflectionCount;

  if (p.useDoubleNucleationSaturation)
  {
    while ((not occupiedNCs.empty()) and occupiedNCs.front() < systemTime + systemTimeOffset - p.occupancytimeNC)
      occupiedNCs.pop_front();
    m.nOccupiedNCs = occupiedNCs.size();
  }
  else
  {
    m.nOccupiedNCs = 0;
  }

	double realVplus;
	double realNuc;
	double freeNucFraction;
  double localCat;
  double majorityRes;
#ifdef BAND_CAT
  localCat = p.catMin;
#else
  localCat = p.kCat;
#endif
	if (p.restrictedPool == 0)
		realVplus = p.vPlus;
	else
		realVplus = p.vPlus*(1.0 - totalLength/(p.poolDensity*geometry->area));
	if (p.ellipseReducedFreeRate) {
		freeNucFraction = p.nucleationHalfIsotropicDensity/(p.nucleationHalfIsotropicDensity + m.lengthDensity);
		realNuc = p.kNuc*(p.ellipseReducedFreeRateAcceptFraction * freeNucFraction + 1.-freeNucFraction);
	}
	else 
    realNuc = p.kNuc;
  switch (p.extraRescueFunction) {
    case r_cos: // 1/2 + 1/2 *cos(2 (theta - kResExtraAngle))
      majorityRes = p.kRes + p.kResExtraMax *( 0.5 + 0.5*cos(2*(m.order.S2angle - p.kResExtraAngle)));
      break;
    default:
      majorityRes = p.kRes;
      break;
  }

#ifdef BAND_CAT
  m.G_effAdjust_band =  (p.kRes/(-p.vMin+p.vTM) - p.kCat[0]/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);
  m.G_effAdjust_gap =  (p.kRes/(-p.vMin+p.vTM) - p.kCat[1]/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);
#endif
  m.G_effAdjust_min =  (p.kRes/(-p.vMin+p.vTM) - localCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);
	m.G_effAdjust_max =  ((p.kRes+p.kResExtraMax)/(-p.vMin+p.vTM) - localCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.); //for BAND_CAT: maybe better change to p.catMax...
	m.G_effAdjust_S2angle =  ((majorityRes)/(-p.vMin+p.vTM) - localCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.); 
	m.G_effAdjust_normal =  (p.kRes/(-p.vMin+p.vTM) - localCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);
	m.G_effAdjust_special =(p.kRes/(-p.vMin+p.vTM) - p.catastropheMultiplier*localCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);
	
	m.G_effMeasured = (m.G_effAdjust_normal)*pow( (p.kNuc)*(m.time)*(geometry->area)/totalNucleationCount, 1./3.);
	//cout << m.time << "	" << totalNucleationCount << "	" << m.G_effAdjust_normal << "	" << totalNucleationCount/((m.time)*(geometry->area)) << endl;

	measurementHistory.push_back(m);

    // The histogram collection would be a good candidate for a separate event (other than measurements)

	if (p.angleHistogramBins)
    {
		h = new double[p.angleHistogramBins];
		geometry->calculateHistogram(h);
		angleHistory.push_back(h);
		if ( p.angleHistogramOptical )
		{
		  h = new double[p.angleHistogramBins];
		  geometry->calculateHistogram(h,true);
		  angleHistoryOptical.push_back(h);
		}
    }


	if (p.hiresLengthHistogramBins)
	{
		// save a microtubule length histogram
		mt = growing_mts.first();
		while (mt != NULL)
		{
			mtLengthHistogram.insert(mt->length());
			mt = mt->next();
		}
		mt = shrinking_mts.first();
		while (mt != NULL)
		{
			mtLengthHistogram.insert(mt->length());
			mt = mt->next();
		}
	}

	if (p.loresLengthHistogramBins && p.loresAngleHistogramBins)
	{
		// save a segment angle/length histogram
		Segment* seg;
		double templength;
		Segment* tempseg;
		mt = growing_mts.first();
		while (mt != NULL)
		{
			seg = mt->segments.first();
			while (seg != NULL)
			{
				if (seg->trajectory->base.region->accountingType == r_accounting_normal)
				{
					if ((seg == mt->segments.first()) 
							|| ((seg->startItr != seg->trajectory->wallBegin()) && (seg->startItr != seg->trajectory->wallEnd())) 
							|| (seg->previous()->trajectory->base.region->accountingType == r_accounting_dontcount))
					{
						templength = seg->length();
						tempseg = seg->next();
						while ((tempseg != NULL) && ((tempseg->startItr == tempseg->trajectory->wallBegin()) || (tempseg->startItr == tempseg->trajectory->wallEnd())))
						{
							templength += tempseg->length();
							tempseg = tempseg->next();
						}
						segAngleLengthHistogram.insert(seg->trajectory->base.angle, templength);
					}
				}
				seg = seg->next();
			}
			mt = mt->next();
		}
		mt = shrinking_mts.first();
		while (mt != NULL)
		{
			seg = mt->segments.first();
			while (seg != NULL)
			{
				if (seg->trajectory->base.region->accountingType == r_accounting_normal)
				{
					if ((seg == mt->segments.first()) 
							|| ((seg->startItr != seg->trajectory->wallBegin()) && (seg->startItr != seg->trajectory->wallEnd())) 
							|| (seg->previous()->trajectory->base.region->accountingType == r_accounting_dontcount))
					{
						templength = seg->length();
						tempseg = seg->next();
						while ((tempseg != NULL) && ((tempseg->startItr == tempseg->trajectory->wallBegin()) || (tempseg->startItr == tempseg->trajectory->wallEnd())))
						{
							templength += tempseg->length();
							tempseg = tempseg->next();
						}
						segAngleLengthHistogram.insert(seg->trajectory->base.angle, templength);
					}
				}
				seg = seg->next();
			}
			mt = mt->next();
		}
	}

	histogramAverageCount++;
	if (histogramAverageCount >= p.histogramAverageSamples)
	{
		if (mtLengthHistogram.valid())
			mtLengthHistogram.save();
		if (mtLifetimeHistogram.valid())
			mtLifetimeHistogram.save();
		if (segAngleLengthHistogram.valid())
			segAngleLengthHistogram.save();
		if (segAngleLifetimeHistogram.valid())
			segAngleLifetimeHistogram.save();
		histogramAverageCount = 0;
	}

	cout << "[t=" << setprecision(9) << m.time << setprecision(6) << "\t\t" << m.trajectories << " trajectories \t" << m.segments << " segments]\n" \
		<< "density [real]: " << m.lengthDensity << "\t\t [optical]: " << m.opticalDensity << "\n" \
		<< "length/MT     : " << m.averageLength << "\tMT#: " << m.numberOfMTs << "\tsegments/MT: " << m.segmentsPerMT << "\n" \
		<< "S2 [real]     : " << m.order.S2 << "\tS2 angle:" << m.order.S2angle<< "\n"\
		<< "S2 [optical]  : " << m.order.S2Opt << "\tS2 angle:" << m.order.S2angleOpt << "\n"\
		<< "S4 [real]     : " << m.order.S4 << "\tS4 angle:" << m.order.S4angle<< "\n"\
		<< "S4 [optical]  : " << m.order.S4Opt << "\tS4 angle:" << m.order.S4angleOpt << "\n"\
		<< "R  [real]     : " << m.order.R << " [" << m.order.Rdirector[0] << ", " << m.order.Rdirector[1] << ", " << m.order.Rdirector[2] << "]\n"\
		<< "G_eff_adjust  : " << m.G_effAdjust_normal << " [on regular surfaces]\n"\
		<< "G_eff_measured: " << m.G_effMeasured << " [on regular surfaces]\n\n";


	if ((p.discreteAngleNumber) && ((geometry->type == g_periodic) || (geometry->type == g_grid)))
	{
		al = new double[p.discreteAngleNumber];
		an = new int[p.discreteAngleNumber];
		switch (geometry -> type)
		{
			case g_periodic:
				static_cast<Periodic*>(geometry)->calculateDiscreteAngleLengths(al,an);
				break;
			case g_grid:
				static_cast<Grid*>(geometry)->calculateDiscreteAngleLengths(al,an);
				break;
		}
		angleLengthHistory.push_back(al);
		angleNumberHistory.push_back(an);

		cout << "angle: length/MT";
		for (i=0; i<p.discreteAngleNumber; i++)
			cout << '\t' << p.nucleationAngles[i] << ": " << al[i];
		cout << "\nangle: #MTs";
		for (i=0; i<p.discreteAngleNumber; i++)
			cout << "\t\t" << p.nucleationAngles[i] << ": " << an[i];
		cout << "\n\n";	
	}	
    if (p.spatialHistogramType == sh_x || p.spatialHistogramType == sh_xy)
        geometry->getOneDmeasurement((*histX),m.time);
    if (p.spatialHistogramType == sh_y || p.spatialHistogramType == sh_xy)
        geometry->getOneDmeasurement((*histY),m.time);

  if (p.outputNucPos == "X" or p.outputNucPos == "XY")
  {
    writeNucleationPositions(p.outputDir, nucleationXpositions, systemTime + systemTimeOffset);
    nucleationXpositions.clear();
  }
  if (p.outputNucPos == "Y" or p.outputNucPos == "XY")
  {
    writeNucleationPositions(p.outputDir, nucleationYpositions, systemTime + systemTimeOffset);
    nucleationYpositions.clear();
  }
        
	return;
}


void System::writeMeasurementsToFile(int numberToKeep)
{
	int i,j, iMax;

	iMax = measurementHistory.size()-numberToKeep;
	for (i=0; i<iMax; i++)
	{
		measurementFile << measurementHistory[0];
    if (p.regionalOutputQuantities and p.geometry == g_gridcylinder)
    {
      writeRegionOrder(measurementHistory[0].local_order, p.outputDir, static_cast<GridCylinder*>(geometry)->number, measurementHistory[0].time);
    }
		measurementHistory.pop_front();
	}


    if (p.angleHistogramBins)
    {
	iMax = angleHistory.size()-numberToKeep;
    	for (i=0; i<iMax; i++)
        {
        	for (j=0; j<p.angleHistogramBins; j++)
            {
            	angleHistogramFile << *(angleHistory[0]+j) << "\t";
            }
            delete angleHistory[0];
            angleHistory.pop_front();
            angleHistogramFile << "\n";
        }
    if ( p.angleHistogramOptical ) {
      iMax = angleHistoryOptical.size()-numberToKeep;
      for (i=0; i<iMax; i++)
      {
        for (j=0; j<p.angleHistogramBins; j++)
        {
          angleHistogramOpticalFile << *(angleHistoryOptical[0]+j) << "\t";
        }
        delete angleHistoryOptical[0];
        angleHistoryOptical.pop_front();
        angleHistogramOpticalFile << "\n";
      }
    }
    }

    if ((p.discreteAngleNumber) && ((geometry->type == g_periodic) || (geometry->type == g_grid)))
    {

	iMax = angleLengthHistory.size()-numberToKeep;
        for (i=0; i<iMax; i++)
        {
            for (j=0; j<p.discreteAngleNumber; j++)
            {
               angleLengthFile << *(angleLengthHistory[0]+j) << "\t";
            }
            delete angleLengthHistory[0];
            angleLengthHistory.pop_front();
            angleLengthFile << "\n";
        }
		
	iMax = angleNumberHistory.size()-numberToKeep;
        for (i=0; i<iMax; i++)
        {
            for (j=0; j<p.discreteAngleNumber; j++)
            {
                angleNumberFile << *(angleNumberHistory[0]+j) << "\t";
            }
            delete angleNumberHistory[0];
            angleNumberHistory.pop_front();
            angleNumberFile << "\n";
         }
    }
    
    #ifdef MTBASED_NUCLEATION_PROBABILITY
    for (int i=0; i<nearbyDensityVector.size(); ++i)
    {
    	MTbasedProbDensityFile << regionalDensityVector[i] << "\t" << nearbyDensityVector[i] << "\t" << globalDensityVector[i] << "\n";
    }
    for (int i=0; i<nearbyDensityVectorUnbound.size(); ++i)
    {
    	unboundProbDensityFile << regionalDensityVectorUnbound[i] << "\t" << nearbyDensityVectorUnbound[i] << "\t" << globalDensityVectorUnbound[i] << "\n";
    }
    #endif
    
    return;
}



ExtensibleHistogram::ExtensibleHistogram(int numberOfBins, double initialRange) :
		histogram(NULL)
{
	initialize(numberOfBins, initialRange);
	return;
}

ExtensibleHistogram::~ExtensibleHistogram()
{
	if (histogram != NULL)
		delete histogram;
}



void ExtensibleHistogram::initialize(int numberOfBins, double initialRange)
{
	range = initialRange;
	totalCount = 0;
	bins = numberOfBins & (-2);	// to make sure it's an even number

	if (histogram != NULL)
		delete histogram;

	if (numberOfBins != 0)
		histogram = new int[bins];
	else
		histogram = NULL;

	for (int i=0 ; i<bins ; i++)
		histogram[i] = 0;
	return;
}

void ExtensibleHistogram::insert(double in)
{
	while (in >= range)
		extendRange();
	totalCount++;
	histogram[static_cast<int>(bins*in/range)]++;
	return;
}
void ExtensibleHistogram::flush(ostream& out)
{
	int i;
	int lastNonZero = 0;
	out << " " << range;
	for (i=0 ; i<bins ; i++)
	{
		out << " " << histogram[i];
		if (histogram[i] != 0)
			lastNonZero = i;
		histogram[i] = 0;
	}
	out << "\n";

	double dx = range/static_cast<double>(bins);
	while (dx*(lastNonZero+1) < 0.5*range)
		range /= 2.0;

	totalCount = 0;
	return;
}

void ExtensibleHistogram::extendRange()
{
	int i;
	for (i=0 ; i< bins/2 ; i++)
		histogram[i] = histogram[2*i] + histogram[2*i+1];
	for (i=bins/2 ; i< bins ; i++)
		histogram[i] = 0;
	range *= 2.0;
	return;
}


SingleHistogram::SingleHistogram(System* s)
		: bins(0),
		  h(0),
		  system(s)
{
	return;	
}

void SingleHistogram::setBins(int b)
{
	if (b != bins)
	{
		if (!h.empty())
			save();
		bins = b;
		h.initialize(b);
	}
	return;
}

void SingleHistogram::setFileName(string f)
{
	if (!h.empty())
		save();

	if (file.is_open())
		file.close();
	file.open(f.c_str());
	if (!file)
	{
		cerr << "Unable to open histogram output file [" << f << "]. Aborting.\n";
		return;
	}
}

void SingleHistogram::closeFile()
{
	if (file.is_open())
		file.close();
	return;
}

SingleHistogram::~SingleHistogram()
{
	if (file.is_open())
		file.close();
	return;
}

void SingleHistogram::insert(double in)
{
	if (bins != 0)
		h.insert(in);
	return;
}

void SingleHistogram::save()
{
	if (file.is_open())
	{
		file << "time " << system->systemTime + system->systemTimeOffset << "\n";
		file << "numberOfAngles 1\n";
		file << "numberOfBins " << bins << "\n";
		file << "0 ";
		h.flush(file);
	}
	else
	{
		cerr << "SingleHistogram::save() : file is closed. Cannot save.\n";
	}
	return;
}

MultiAngleHistogram::MultiAngleHistogram(System* s) :
			angles(0),
			bins(0),
			h(NULL),
			system(s)
{
	return;
}

void MultiAngleHistogram::setBins(int numAngles, int b)
{
	if (((angles != numAngles) || (bins != b)) && (numAngles != 0) && (b != 0))
	{
		for (int i=0 ; i<angles ; i++)
		{
			if (!h[i]->empty())
			{
				save();
				break;
			}
		}
		angles = numAngles;
		bins = b;

		if (h != NULL)
		{
			for (int i=0 ; i<angles ; i++)
				delete h[i];
			delete h;
		}
		h = new ExtensibleHistogram *[angles];
		for (int i=0 ; i<angles ; i++)
			h[i] = new ExtensibleHistogram(b);
	}
	return;
}

void MultiAngleHistogram::setFileName(string f)
{
	for (int i=0 ; i<angles ; i++)
	{
		if (!h[i]->empty())
		{
			save();
			break;
		}
	}
	if (file.is_open())
		file.close();
	file.open(f.c_str());
	if (!file)
	{
		cerr << "Unable to open histogram output file [" << f << "]. Aborting.\n";
		exit(-1);
	}
	return;
}

void MultiAngleHistogram::closeFile()
{
	if (file.is_open())
		file.close();
	return;
}


MultiAngleHistogram::~MultiAngleHistogram()
{
	if (file.is_open())
		file.close();
	return;
}

void MultiAngleHistogram::insert(double angle, double length)
{
	if ((bins != 0) && (angles !=0))
	{
		int bin = static_cast<int>(((angle/PI + 2.0)*angles + 0.5))%angles;
		h[bin]->insert(length);
	}
	return;
}

void MultiAngleHistogram::save()
{
	if (file.is_open())
	{
		file << "time " << system->systemTime + system->systemTimeOffset << "\n";
		file << "numberOfAngles " << angles << "\n";
		file << "numberOfBins " << bins << "\n";
		for (int i=0 ; i<angles ; i++)
		{
			file << i*PI/angles << " ";
			h[i]->flush(file);
		}
	}
	else
	{
		cerr << "MultiAngleHistogram::save() : file is closed. Cannot save.\n";
	}
	return;
}

OneDSpatialMeasurement::OneDSpatialMeasurement(int n, bool caps, OrientationType oo, double len, double ts, bool wr) : bins(n),includeCaps(caps),ori(oo),binSurface(ts/n),writeRaw(wr){
    double dlen;
    dlen = len/n;
    cuts = new double[n];
    cuts[n-1] = VERY_LARGE;
    for (int i=0 ; i < n-1; i++)
    {
        cuts[i] = -0.5*len + dlen*(i+1)  ;
    }
    if (includeCaps){
        opHist = new OrderParametersRaw[n+2];
        nHist = new double[n+2];
        cuts[n] = 0.;
        cuts[n+1]= len;
    }
    else 
    {
        opHist = new OrderParametersRaw[n];
        nHist = new double[n];
    }
    //for (int i = 0 ; i< n + 2*includeCaps ; i++) {
        //cout << cuts[i] << "\t";
    //}
    //cout << "\n";
    return;
}

OneDSpatialMeasurement::~OneDSpatialMeasurement() {
	delete[] cuts;
	delete[] opHist;
	delete[] nHist;
	return;
}

void OneDSpatialMeasurement::writeHeader() {
    if (writeRaw){
        file << "#position\tMTcount\tdensity\tS2\tS2_angle\tS4\tS4_angle\tdensityOpt\tS2Opt\tS2Opt_angle\tS4Opt\tS4Opt_angle\tR\tR_x\tR_y\tR_z";
        file << "\tsi2\tsi4\tco2\tco4\tlocalL\tsi2Opt\tsi4Opt\tco2Opt\tco4Opt\tlocalLOpt\tQxx\tQxy\tQxz\tQyy\tQyz\tQzz\tisoWeights0\tisoWeights1\tisoWeights2\n";
    }
    else
        file << "#position\tMTcount\tdensity\tS2\tS2_angle\tS4\tS4_angle\tdensityOpt\tS2Opt\tS2Opt_angle\tS4Opt\tS4Opt_angle\tR\tR_x\tR_y\tR_z\n";
}

void OneDSpatialMeasurement::writeLine(double pos, int i) {
    double R, Rdirector[3];
    R = opHist[i].extractR(Rdirector);
    file << pos << "\t" << nHist[i] << "\t" \
        << opHist[i].localL/binSurface << "\t" \
        << sqrt(opHist[i].si2*opHist[i].si2 + opHist[i].co2*opHist[i].co2) << "\t" \
        << 0.5*atan2(opHist[i].si2,opHist[i].co2) << "\t" \
        << sqrt(opHist[i].si4*opHist[i].si4 + opHist[i].co4*opHist[i].co4) << "\t" \
        << 0.25*atan2(opHist[i].si4,opHist[i].co4) << "\t" \
        << opHist[i].localLOpt/binSurface << "\t" \
        << sqrt(opHist[i].si2Opt*opHist[i].si2Opt + opHist[i].co2Opt*opHist[i].co2Opt) << "\t" \
        << 0.5*atan2(opHist[i].si2Opt,opHist[i].co2Opt) << "\t" \
        << sqrt(opHist[i].si4Opt*opHist[i].si4Opt + opHist[i].co4Opt*opHist[i].co4Opt) << "\t" \
        << 0.25*atan2(opHist[i].si4Opt,opHist[i].co4Opt) << "\t" \
        << R << "\t" << Rdirector[0] << "\t" << Rdirector[1] << "\t" << Rdirector[2]; 
    if (writeRaw) {
        file << "\t" << opHist[i].si2 * opHist[i].localL << "\t" \
        << opHist[i].si4 * opHist[i].localL << "\t" \
        << opHist[i].co2 * opHist[i].localL << "\t" \
        << opHist[i].co4 * opHist[i].localL << "\t"  << opHist[i].localL << "\t"  \
        << opHist[i].si2Opt * opHist[i].localLOpt << "\t" \
        << opHist[i].si4Opt * opHist[i].localLOpt << "\t" \
        << opHist[i].co2Opt * opHist[i].localLOpt << "\t" \
        << opHist[i].co4Opt * opHist[i].localLOpt << "\t"  << opHist[i].localLOpt << "\t"  \
        << opHist[i].Qxx * opHist[i].localL << "\t" \
        << opHist[i].Qxy * opHist[i].localL << "\t" \
        << opHist[i].Qxz * opHist[i].localL << "\t" \
        << opHist[i].Qyy * opHist[i].localL << "\t" \
        << opHist[i].Qyz * opHist[i].localL << "\t" \
        << opHist[i].Qzz * opHist[i].localL << "\t" \
        << opHist[i].isoWeights[0] * binSurface << "\t" \
        << opHist[i].isoWeights[1] * binSurface << "\t" \
        << opHist[i].isoWeights[2] * binSurface;
    }
    file << "\n";
    return;
}

void OneDSpatialMeasurement::writeToFile(double t) {
    // IMPORTANT! assumes that ::process() is called just before.
    double pos;
    //file << "#Time "<<setprecision(9) << t << " s\n"; 
    file << "#Time " << t << " s\n"; 
   if (includeCaps ) 
       writeLine(0.,bins);
   for (int i=0; i<bins ; i++) {
       if (i == 0)
           pos=0.5*(cuts[0]+cuts[1])-(cuts[1]-cuts[0]);
       else if (i == bins -1)
           pos=0.5*(cuts[bins-2]+cuts[bins-3]) + cuts[1] - cuts[0]; // assume equally spaced cuts!
       else
           pos=0.5*(cuts[i]+cuts[i-1]);
       writeLine(pos,i);
   } 
   if (includeCaps ) 
       writeLine(cuts[bins+1],bins+1);
   file << "\n";
   return;
}

void OneDSpatialMeasurement::process() {
    // IMPORTANT! Only call before writing to file. Needs reset before meaningful new measurement can be made.
    for (int i= 0 ; i < bins + 2*includeCaps; i++ ){
        if (opHist[i].localL > ZERO_CUTOFF)
        {
            opHist[i].si2 /= opHist[i].localL;
            opHist[i].si4 /= opHist[i].localL;
            opHist[i].co2 /= opHist[i].localL;
            opHist[i].co4 /= opHist[i].localL;
            opHist[i].si2Opt /= opHist[i].localLOpt;
            opHist[i].si4Opt /= opHist[i].localLOpt;
            opHist[i].co2Opt /= opHist[i].localLOpt;
            opHist[i].co4Opt /= opHist[i].localLOpt;

            opHist[i].Qxx /= opHist[i].localL;
            opHist[i].Qxy /= opHist[i].localL;
            opHist[i].Qxz /= opHist[i].localL;
            opHist[i].Qyy /= opHist[i].localL;
            opHist[i].Qyz /= opHist[i].localL;
            opHist[i].Qzz /= opHist[i].localL;
        }
            opHist[i].isoWeights[0] /= binSurface;
            opHist[i].isoWeights[1] /= binSurface;
            opHist[i].isoWeights[2] /= binSurface;
    }
    return;
}
void OneDSpatialMeasurement::reset(void) 
{
    for (int i=0; i<bins+2*includeCaps; i++) 
    {   
        nHist[i] = 0.;
        opHist[i].si2 =0.;
        opHist[i].si4 =0.;
        opHist[i].co2 =0.;
        opHist[i].co4 =0.;
        opHist[i].localL =0.;
        opHist[i].si2Opt =0.;
        opHist[i].si4Opt =0.;
        opHist[i].co2Opt =0.;
        opHist[i].co4Opt =0.;
        opHist[i].localLOpt =0.;
        opHist[i].Qxx =0.;
        opHist[i].Qxy =0.;
        opHist[i].Qxz =0.;
        opHist[i].Qyy =0.;
        opHist[i].Qyz =0.;
        opHist[i].Qzz =0.;
        opHist[i].isoWeights[0] =0.;
        opHist[i].isoWeights[1] =0.;
        opHist[i].isoWeights[2] =0.;
      } 
    return;
}

void OneDSpatialMeasurement::setFileName(string f)
{
	if (file.is_open())
		file.close();
	file.open(f.c_str());
	if (!file)
	{
		cerr << "Unable to open histogram output file [" << f << "]. Aborting.\n";
		return;
	}
}

void OneDSpatialMeasurement::closeFile()
{
	if (file.is_open())
		file.close();
	return;
}

