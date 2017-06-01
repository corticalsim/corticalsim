/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
				<< m.G_effAdjust_normal << "\t" << m.G_effAdjust_special << "\r\n";
}

void writeMeasurementDescriptors(ostream& o)
{
	o << "time\tdensity"\
		<< "\t<l>" \
		<< "\tS2\tS2 angle\tS4\tS4 angle\t#growing\t#shrinking\t#segments\t#trajectories"\
		<< "\tzippering events\tcrossover events\tinduced catastrophe events\t"\
		<< "valid deterministic events\tinvalid deterministic events\tstochastic events"\
		<< "\toptical density\tmicrotubules\tsegments per MT\trandom severing events\tintersection severing events"\
		<< "\tS2 opt\tS2 opt angle\tS4 opt\tS4 opt angle\toccupied intersections"\
		<< "\tR\tR_x\tR_y\tR_z"\
		<< "\tG_eff_adjusted normal\tG_eff_adjusted special\r\n";
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

	mtLengthHistogram.closeFile();
	mtLifetimeHistogram.closeFile();
	segAngleLengthHistogram.closeFile();
	segAngleLifetimeHistogram.closeFile();

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


	m.zipperCount = totalZipperCount;
	m.crossoverCount = totalCrossoverCount;
	m.inducedCatastropheCount = totalInducedCatastropheCount;
	m.validDEventCount = totalValidDEventCount;
	m.invalidDEventCount = totalInvalidDEventCount;
	m.sEventCount = totalSEventCount;
	m.lengthSeveringCount = totalLengthSeveringCount;
	m.intersectionSeveringCount = totalIntersectionSeveringCount;
	m.occupiedIntersectionCount = OccupiedIntersectionList.size();

	double realVplus;
	double realNuc;
	double freeNucFraction;
  realVplus = p.vPlus;
	if (p.ellipseReducedFreeRate) {
		freeNucFraction = p.nucleationHalfIsotropicDensity/(p.nucleationHalfIsotropicDensity + m.lengthDensity);
		realNuc = p.kNuc*(p.ellipseReducedFreeRateAcceptFraction * freeNucFraction + 1.-freeNucFraction);
	}
	else 
		realNuc = p.kNuc;
	m.G_effAdjust_normal =  (p.kRes/(-p.vMin+p.vTM) - p.kCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);
	m.G_effAdjust_special =(p.kRes/(-p.vMin+p.vTM) - p.catastropheMultiplier*p.kCat/(realVplus - p.vTM))*pow((2.0*(-p.vMin+p.vTM)*(realVplus-p.vTM)*(realVplus - p.vTM))/(realNuc*(-p.vMin+realVplus)*realVplus), 1./3.);

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
		<< "G_eff_adjust  : " << m.G_effAdjust_normal << " [on regular surfaces]\n\n";


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

	return;
}


void System::writeMeasurementsToFile(int numberToKeep)
{
	int i,j, iMax;

	iMax = measurementHistory.size()-numberToKeep;
	for (i=0; i<iMax; i++)
	{
		measurementFile << measurementHistory[0];
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


