#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <cmath>
#include <cstdlib>
#include <limits>



#pragma warning(disable:981)
//#pragma warning(disable:177)
#pragma warning(disable:522)
#pragma warning(disable:383)


#include "corticalSim.h"


string GeometryTypeText[] = {"periodic", "grid", "wormhole", "box", "cylinder", "gridcylinder", "spherocylinder", "sphere", "pancake"};
string RegionTypeText[] = {"rectangle", "disc", "dome"};
string OrientationTypeText[] = {"x", "y", "z"};
string SpatialHistogramTypeText[] = {"none", "x", "y", "xy" };
string NucleationTypeText[] = {"isotropic", "aster", "biased", "discreteAngles", "ellipse" };
string NucleationBiasTypeText[] = {"cross", "ellipsePolar", "ellipseApolar" };
//string NucleationTypeText[] = {"isotropic", "aster", "biased", "discreteAngles", "chanLloyd", "chanLloydRandomPosition", "chanLloydIsotropic", "ellipse" };
string InteractionTypeText[] = {"zipFirst", "catFirst", "minimalFourier"};
string BundleTypeText[] = {"simple", "sticky", "noZip", "multiCollision", "Ncollision" };
string PatternTypeText[] = {"bar", "cross", "band", "singleBand"};
string extraRescueTypeText[] = {"cos"}; // for future: option to change function for extra rescue


Parameters::Parameters(System* s) :
		system(s),

		vPlus(0.05),
		vMin(-0.08),
		vTM(0.01),
		kSev(0.001),
		kCross(0.001),
		kResExtraMax (0.),
		kResExtraAngle (0.),                  // angle of highest rescue rate
		extraRescueFunction (r_cos), 
		
		// parameters with Aster creation
		numberOfRays(3),
		maxMetaTrajectoryLength(1.001),
		diffusionCoefficient(0.013),
		rvRejectionUnbound(0.),
		rvRejectionBound(0.),
		unboundNucleationRate(0.002),
		
#ifdef BAND_CAT
    kCatOrient (o_x),  // Orientation: along x, y, or z axis (z not implemented!)
    nSpirals(0),
    wrapLength(1.),
    projectedBand(0),
    projectedPeriod(0),
    calcNucBiasAngle(0),
    catMin(0.), // theoretical checks (4x)
    catMax(0.),
    Gmin(0.),
    Gmax(0.),
    catPattern(p_band), // default: bands or spirals (both p_band)
    redistributeNucleationOverBands(false),
    reducedGapNucleation(false),
    gapNucleationAcceptFraction(1.),
#else
		kCat(0.005),
		Gmax(0.),                   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
#ifdef VAR_CAT
    kCatOrient (o_y),  // Orientation: along x, y, or z axis (z not implemented!)
    kCatAlpha (0.), // 3 parameters used for sin-gradient of kCat rates. 
    kCatBands (0),  // n bands over total y length of geometry
    invKCatPeriod (0.), // bookkeeping: 2 * pi * yLength / kCatBands
    catMin(0.), // theoretical checks (4x)
    catMax(0.),
    Gmin(0.),
    Gmax(0.),
    //redistributeNucleationOverBands(false),  // RECENT EDIT: was in wrong ifdef
#endif
		kRes(0.007),
		poolDensity(10.0),
		treadmillingEnabled(1),
		severingEnabled(0),
		crossSeveringEnabled(0),
		catastropheMultiplier(1.0),
		forbiddenZones(0),

		edgeCatastropheEnabled(0),
		edgeCatastropheSmooth(0),
		pCatRegularEdge(0.0),
		pCatSpecialEdge(0.0),

		kNuc(0.0045),
		nucleationType(nuc_aster),
		nucleationBiasType(nbias_cross),  // default to maintain previous meaning of nuc_biased; only has a meaning with biased or ellipse nucleation
		nucleationAlpha(0.0),
		nucleationBiasAngle(0.0),
		nucleationAngles(NULL),

		interactionType(int_zipFirst),
		zipperingEnabled(1),
		catastrophesEnabled(1),
		proportionalCatastrophes(0),
		inducedCatastropheFraction(0.5),
		catStartAngle(0),
		zipFraction(1.0),
		magicAngle(40),
		c0Value(0.75),
		z0Value(0.5),

		crossSeveringTop(1),
		crossSeveringTopFraction(0.5),
		crossSeveringStartAngle(0),

		restrictedPool(0),
		
		bundleType(bdl_simple),
		discreteAngleNumber(0),

		
		geometry(g_periodic),
		geomParam1(35.4),
		geomParam2(70.8),
		geomParam3(50),
		
		//seed(time(0)), This code is intuitive, but not guaranteed to work. See below for new initialisation code.
		stopTime(36000),
		measurementInterval(600),
		densityLimit(VERY_LARGE),
		wallClockLimit(VERY_LARGE),
		memoryLimit(VERY_LARGE),
		
		outputDir("."),
		createSubdir(1),
		newParameterReadInterval(VERY_LARGE),
		newParameterFile(""),


		movieEnabled(0),			//make sure movie file is always created
		movieFrameInterval(-1),
		angleHistogramBins(0),			
		angleHistogramOptical(0),			
		loresLengthHistogramBins(0),
		loresAngleHistogramBins(0),
		loresLifetimeHistogramBins(0),
		hiresLengthHistogramBins(0),
		hiresLifetimeHistogramBins(0),
		histogramAverageSamples(1),
        spatialHistogramType(sh_none),
        spatialHistogramBinsX(0),
        spatialHistogramBinsY(0),
        spatialHistogramCountCaps(false),
        spatialHistogramWriteRaw(false),
        regionalOutputQuantities(false),

		c0calc(0),
		x0calc(0),
		z0calc(0),
		G_0(0),
		adjustedG(0),
		lengthFactor(0),
    //preSeededNucleation(false),           // a fixed number of seeds present at the beginning of the simulation, they fire once with a certain rate and can have different properties than other (free!!) nucleations. Always unbound. // Not a separate bool needed. 
    preSeededSeedDensity(0.),             // sets number of seeds (Density * Area)
    preSeededRate(0.),                    // fire rate of seeds (a seed can start growing only once)
    preSeededType(nuc_isotropic),         // options available: isotropic and biased. 
    preSeededAlpha(0.),                   // for biased seeded nucleations.
    nucleationHalfIsotropicDensity(1.),		// When less segments are present, use random isotropic nucleation.
    ellipseReducedFreeRate(0),                  // When true, a fraction of free nucleation events is rejected.
    ellipseReducedFreeRateAcceptFraction(1.),       // This fraction of free nucleation events is accepted. It is impossible to increase the free nucleation rate by choosing this >1.
    ellipseEpsilon(0.),				// eccentricity of ellipse: [0,1) (0 = isotropic, 1 = delta function (putting 1 will lead to errors!)
    ellipseEpsilonAlongMT(-1.),// like epsilon, but specific for parallel and antiparallel nucleations (default -1 will be set to the value in ellipseEpsilon)
    ellipseForwardAlongMT(0),			// forward/backward nucleation not elliptic, but along MT
		ellipseLeftFraction(0.),			// fraction of nucleations rotated left over ellipseSidewaysAngle degrees
		ellipseRightFraction(0.),			// fraction of nucleations rotated right over ellipseSidewaysAngle degrees
		ellipseBackwardFraction(0.),			// fraction of nucleations rotated 180 degrees
		ellipseSidewaysAngle(40.),			// mean angle of sideways nucleations (+/- ellipseSidewaysAngle degrees)
    useRegionalNucleationSaturation(false),   // apply saturation of bound nucleations regionally, with a fixed unbound nucleation rate
    f_unbound(0.),     // fraction of kNuc that will be unbound if useRegionalNucleationSaturation == true
    useDoubleNucleationSaturation(false),   // apply saturation of bound nucleations regionally, with max and base rates depending on global NC availability
    rn_targetratio(0.),     // ratio between rnmax at the start and rnmax for an established homogeneous array
    occupancytimeNC(0.),     // time a nucleation complex spends occupied at the membrane upon nucleating before being available for a new nucleation
    nNCmax(0),     // extra par, not read but derived from others
    rn_base0(0.),     // extra par, not read but derived from others
    gridAspectratio(1.),     // aspect ratio of grid cells for gridcylinder or grid geometry
    useAbsoluteBidirectionalNucBias(false),   // turn on to give unbound nucleations a specific orientation
    absoluteNucBiasAngle(0.),                  // angle of the orientation that unbound nucleations will get
    nucBiasVariance(0.),                  // angle of the orientation that unbound nucleations will get
    useMTdeflection(true),   // have MTs change direction every deflectionStepsize by a uniform angle between + and - deflectionMaxAngle
    deflectionStepsize(3.),     // length at which next deflection is to occur
    deflectionMaxAngle(PI/9.),     // deflection angles are drawn uniformly from [-deflectionAngle, +deflectionAngle]
    deflectionMinAngle(PI/180.),     // deflection angles are drawn uniformly from [-deflectionAngle, +deflectionAngle] excluding [-minAngle,minAngle]
    deflectionBundleCompensation(true),     // reduce deflection rate of MTs in bundles proportional to the number of MTs in the bundle
    deflectionBundleProtection(false),     // no deflection of MTs in bundles
    biasDeflectionRight(0.5),
    useBundleTracking(false),   // MTs in bundle will pick random MT from bundle to follow at intersection
    bundleTrackMaxAngle(PI),   // Maximum angle along which to track MTs in bundle
    outputNucPos("None")   // Specify whether to store nucleation x and y positions
{
// WARNING: Do not use the system pointer inside this function, as it has not been initialized

	// Code for initialisation of 'seed' using current time. Time is given as a time_t type, which is not necessarily
	// expressed in seconds. This code explicitly transforms the current time_t to seconds from 1 Jan 1970. 
	// Based on code from http://www.cplusplus.com/reference/ctime/time/?kw=time 
	struct tm y1970;		// make Unix start date
	y1970.tm_hour = 0;   y1970.tm_min = 0; y1970.tm_sec = 0;
	y1970.tm_year = 70; y1970.tm_mon = 0; y1970.tm_mday = 1;
	time_t timer = time(NULL);
	seed = static_cast<unsigned long>(difftime(timer,mktime(&y1970)));
#ifdef BAND_CAT
		kCat[0] = 0.005;
    kCat[1] = 0.006;
    bandGapWidth[0] = 0.;
    bandGapWidth[1] = 0.;
#endif
	
	return;
}

void Parameters::initialize(const char* pf)
{
	if (!readFromFile(pf, true))
	{
		cout << "Failed to load parameters from file [" << pf << "]. Aborting\n";
		exit(-1);
	}
	verifyParameters();

	calcTheoryParameters();

	return;
}

bool Parameters::reinitialize(const char* pf)
{
	if (readFromFile(pf, false))
	{
		verifyParameters();


		calcTheoryParameters(); 
		return true;
	}
	else
		return false;
}


bool Parameters::writeToFile()
{

	ofstream& of = system->parameterFile;

	of << "# CorticalSim parameters\r\n";

	of << "corticalSim\t" << PROGRAM_VERSION << "\r\n";
	of << "outputDir\t" << outputDir << "\r\n";

	of << "newParameterReadInterval\t" << newParameterReadInterval << "\r\n";
	if (newParameterReadInterval > 0)
		of << "newParameterFile\t" << newParameterFile << "\r\n";


	of << "vPlus\t" << vPlus << "\r\n";
	of << "vMin\t" << vMin << "\r\n";
	of << "vTM\t" << vTM << "\r\n";
	of << "kSev\t" << kSev << "\r\n";
	of << "kCross\t" << kCross << "\r\n";
	
	of << "numberOfRays\t" << numberOfRays << "\r\n";
	of << "maxMetaTrajectoryLength\t" << maxMetaTrajectoryLength << "\r\n";
	of << "diffusionCoefficient\t" << diffusionCoefficient << "\r\n";
	of << "rvRejectionUnbound\t" << rvRejectionUnbound << "\r\n";
	of << "rvRejectionBound\t" << rvRejectionBound << "\r\n";
	of << "unboundNucleationRate\t" << unboundNucleationRate << "\r\n";
	
#ifdef BAND_CAT
	of << "kCat\t" << kCat[0] << "\t" << kCat[1] << "\r\n";
	of << "kCatOrient\t" << OrientationTypeText[kCatOrient] << "\r\n";
  of << "bandGapWidth\t" << bandGapWidth[0] << "\t" << bandGapWidth[1] << "\r\n";
  of << "nSpirals\t" << nSpirals << "\r\n";
	of << "catPattern\t" << PatternTypeText[catPattern] << "\r\n";
  of << "spiralPitch\t" << spiralPitch << "\r\n";
  of << "wrapLength\t" << wrapLength << "\r\n";
  of << "projectedPeriod\t" << projectedPeriod << "\r\n";
  of << "projectedBand\t" << projectedBand << "\r\n";
  of << "calcNucBiasAngle\t" << calcNucBiasAngle << "\r\n";
	of << "redistributeNucleationOverBands\t" << redistributeNucleationOverBands << "\r\n";
  of << "redistributeShift\t" << redistributeShift << "\r\n";
  of << "redistributeShiftNumber\t" << redistributeShiftNumber << "\r\n";
  of << "reducedGapNucleation\t" << reducedGapNucleation << "\r\n";
  of << "gapNucleationAcceptFraction\t" << gapNucleationAcceptFraction << "\r\n";
#else
	of << "kCat\t" << kCat << "\r\n";
#endif

#ifdef VAR_CAT
	of << "kCatAlpha\t" << kCatAlpha << "\r\n";
	of << "kCatBands\t" << kCatBands << "\r\n";
	of << "kCatOrient\t" << OrientationTypeText[kCatOrient] << "\r\n";
#endif
	of << "kRes\t" << kRes << "\r\n";
	of << "kResExtraMax\t" << kResExtraMax << "\r\n";
	of << "kResExtraAngle\t" << kResExtraAngle << "\r\n";
	of << "extraRescueFunction\t" << extraRescueTypeText[extraRescueFunction] << "\r\n";
	of << "poolDensity\t" << poolDensity << "\r\n";
	of << "treadmilling\t" << treadmillingEnabled <<"\r\n";
	of << "severing\t" << severingEnabled <<"\r\n";
	of << "crossSevering\t" << crossSeveringEnabled <<"\r\n";
	of << "catastropheMultiplier\t" << catastropheMultiplier << "\r\n";
	of << "forbiddenZones\t" << forbiddenZones << "\r\n";
	of << "restrictedPool\t" << restrictedPool << "\r\n";

	of << "edgeCatastropheEnabled\t" << edgeCatastropheEnabled << "\r\n";
	of << "edgeCatastropheSmooth\t" << edgeCatastropheEnabled << "\r\n";
	of << "pCatRegularEdge\t" << pCatRegularEdge << "\r\n";
	of << "pCatSpecialEdge\t" << pCatSpecialEdge << "\r\n";


	of << "kNuc\t" << kNuc << "\r\n";
	of << "nucleationType\t" << NucleationTypeText[nucleationType] << "\r\n";
	of << "nucleationBiasType\t" << NucleationBiasTypeText[nucleationBiasType] << "\r\n";
	of << "nucleationAlpha\t" << nucleationAlpha << "\r\n";
	of << "nucleationBiasAngle\t" << nucleationBiasAngle << "\r\n";


	of << "zippering\t" << zipperingEnabled <<"\r\n";
	of << "ind_cat\t" << catastrophesEnabled <<"\r\n";
	of << "proportionalCatastrophes\t" << proportionalCatastrophes << "\r\n";
	of << "induced_cat_fraction\t" << inducedCatastropheFraction << "\r\n";
	of << "cat_start_angle\t" << catStartAngle << "\r\n";
	of << "zipFraction\t" << zipFraction << "\r\n";
	of << "magic_angle\t" << magicAngle << "\r\n";
	of << "c0Value\t" << c0Value << "\r\n";
	of << "z0Value\t" << z0Value << "\r\n";

	of << "interactionType\t" << InteractionTypeText[interactionType] << "\r\n";
	of << "bundleType\t" << BundleTypeText[bundleType] << "\r\n";

	of << "discreteAngleNumber\t" << discreteAngleNumber << "\r\n";	



	of << "random_seed\t" << seed << "\r\n";
	of << "stopTime\t" << stopTime << "\r\n";
	of << "measurementInterval\t" << measurementInterval << "\r\n";
	of << "densityLimit\t" << densityLimit << "\r\n";
	of << "wallClockLimit\t" << wallClockLimit << "\r\n";
	of << "memoryLimit\t" << memoryLimit << "\r\n";
	of << "movieEnabled\t" << movieEnabled << "\r\n";
	if (movieEnabled)
		of << "movieFrameInterval\t" << movieFrameInterval << "\r\n";
	of << "angleHistogramBins\t" << angleHistogramBins << "\r\n";
        of << "angleHistogramOptical\t" << angleHistogramOptical << "\r\n";
		
	of << "hiresLengthHistogramBins\t" << hiresLengthHistogramBins << "\r\n";
	of << "loresLengthHistogramBins\t" << loresLengthHistogramBins << "\r\n";
	of << "loresAngleHistogramBins\t" << loresAngleHistogramBins << "\r\n";
	of << "loresLifetimeHistogramBins\t" << loresLifetimeHistogramBins << "\r\n";
	of << "hiresLifetimeHistogramBins\t" << hiresLifetimeHistogramBins << "\r\n";
	of << "histogramAverageSamples\t" << histogramAverageSamples << "\r\n";
	of << "spatialHistogramType\t" << SpatialHistogramTypeText[spatialHistogramType] << "\r\n";
	of << "spatialHistogramBinsX\t" << spatialHistogramBinsX << "\r\n";
	of << "spatialHistogramBinsY\t" << spatialHistogramBinsY << "\r\n";
	of << "spatialHistogramCountCaps\t" << spatialHistogramCountCaps << "\r\n";
	of << "spatialHistogramWriteRaw\t" << spatialHistogramWriteRaw << "\r\n";
	of << "regionalOutputQuantities\t" << regionalOutputQuantities << "\r\n";


	of << "geometry\t" << GeometryTypeText[geometry] << " ";
	of << geomParam1 << " " << geomParam2 << " " << geomParam3 << "\r\n";

	if (crossSeveringEnabled == 1)
	{
		of << "crossSeveringTop\t" << crossSeveringTop << "\r\n";
		of << "crossSeveringTopFraction\t" << crossSeveringTopFraction << "\r\n";
		of << "crossSeveringStartAngle\t" << crossSeveringStartAngle << "\r\n";
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// add all others things related to nuc_aster
	if (nucleationType == nuc_ellipse)
	{
    of << "ellipseReducedFreeRate\t" <<     ellipseReducedFreeRate << "\r\n";
    of << "ellipseReducedFreeRateAcceptFraction\t" <<       ellipseReducedFreeRateAcceptFraction << "\r\n";
		of << "ellipseEpsilon\t" << 	ellipseEpsilon << "\r\n";
		of << "ellipseEpsilonAlongMT\t" << 	ellipseEpsilonAlongMT << "\r\n";
		of << "ellipseForwardAlongMT\t" << 	ellipseForwardAlongMT << "\r\n";
		of << "ellipseLeftFraction\t" << 	ellipseLeftFraction << "\r\n";
		of << "ellipseRightFraction\t" << 	ellipseRightFraction << "\r\n";
		of << "ellipseBackwardFraction\t" << 	ellipseBackwardFraction << "\r\n";
		of << "ellipseSidewaysAngle\t" << 	ellipseSidewaysAngle << "\r\n";
		of << "nucleationHalfIsotropicDensity\t" << 	nucleationHalfIsotropicDensity << "\r\n";
	}

  if (preSeededRate > ZERO_CUTOFF) {
    of << "preSeededType\t" << NucleationTypeText[preSeededType] << "\r\n";
    of << "preSeededAlpha\t" << preSeededAlpha << "\r\n";
    of << "preSeededRate\t" << preSeededRate << "\r\n";
    of << "preSeededSeedDensity\t" << preSeededSeedDensity << "\r\n";
  }

  if (useRegionalNucleationSaturation){
    of << "useRegionalNucleationSaturation\t" << useRegionalNucleationSaturation << "\r\n";
    of << "f_unbound\t" << f_unbound << "\r\n";
  }

  if (useDoubleNucleationSaturation){
    of << "useDoubleNucleationSaturation\t" << useDoubleNucleationSaturation << "\r\n";
    of << "rn_targetratio\t" << rn_targetratio << "\r\n";
    of << "occupancytimeNC\t" << occupancytimeNC << "\r\n";
    of << "nNCmax\t" << nNCmax << "\r\n";
    of << "rn_base0\t" << rn_base0 << "\r\n";
  }

  if (useMTdeflection){
    of << "useMTdeflection\t" << useMTdeflection << "\r\n";
    of << "deflectionStepsize\t" << deflectionStepsize << "\r\n";
    of << "deflectionMaxAngle\t" << deflectionMaxAngle << "\r\n";
    of << "deflectionMinAngle\t" << deflectionMinAngle << "\r\n";
    of << "deflectionBundleCompensation\t" << deflectionBundleCompensation << "\r\n";
    of << "deflectionBundleProtection\t" << deflectionBundleProtection << "\r\n";
    of << "biasDeflectionRight\t" << biasDeflectionRight << "\r\n";
  }

  of << "useBundleTracking\t" << useBundleTracking << "\r\n";
  of << "bundleTrackMaxAngle\t" << bundleTrackMaxAngle << "\r\n";

  of << "outputNucPos\t" << outputNucPos << "\r\n";

  if (geometry == g_gridcylinder or geometry == g_grid){
    of << "gridAspectratio\t" << gridAspectratio << "\r\n";
  }

  if (useAbsoluteBidirectionalNucBias){
    of << "useAbsoluteBidirectionalNucBias\t" << useAbsoluteBidirectionalNucBias << "\r\n";
    of << "absoluteNucBiasAngle\t" << absoluteNucBiasAngle << "\r\n";
    of << "nucBiasVariance\t" << nucBiasVariance << "\r\n";
  }

	of << "\n# Calculated theory parameters \r\n";
	of << "c0\t" << c0calc << "\r\n";
	of << "z0\t" << z0calc << "\r\n";
	of << "x0\t" << x0calc << "\t(should be: " << 4./PI - z0calc - c0calc << ")\r\n";
	of << "G_0\t" << G_0 << "\r\n";
	of << "adjustedG\t" << adjustedG << "\r\n";
	of << "l0\t" << lengthFactor << "\r\n";
#if defined (VAR_CAT) || defined (BAND_CAT)
	of << "Gmax\t" << Gmax << "\r\n";
	of << "kCatMin\t" << catMin << "\r\n";
	of << "Gmin\t" << Gmin << "\r\n";
	of << "kCatMax\t" << catMax << "\r\n";
#endif
	return true;
}
	
	
template<class T> inline void loadParam(bool& recognized, ifstream& file, string& id, const char* tag, T& target)
{
	if (!recognized)
	{
		if (id == tag)
		{
			file >> target;
			cout << "Set " << tag << " : " << target << "\n";
			recognized = true;
		}
	}
	return;
}	

template<class T> inline void loadParamRandom(bool& recognized, ifstream& file, string& id, const char* tag, T& target)
{
	if (!recognized)
	{
		if (id == tag)
		{
			if (id == "random_seed")
			{
				string temp;
        int maxlen;
				file >> temp;
        maxlen = boost::lexical_cast<string>(numeric_limits<unsigned long long int>::max()).length();
        while (temp.length() > maxlen) 
        {
					temp = temp.substr(1);
        }
        if (temp.length() == maxlen && temp > boost::lexical_cast<string>(ULONG_MAX)) 
        {
					temp = temp.substr(1);
        }
				while (boost::lexical_cast<unsigned long long int>(temp) > numeric_limits<T>::max() ) 
				{
					temp = temp.substr(1);
				}
				target = boost::lexical_cast<T>(temp);
			}
			else
			{
				cout << "Using loadParamRandom for the wrong paramter: " << id << ", quitting.\n";
				exit(-21254);
			}
			cout << "Set " << tag << " : " << target << "\n";
			recognized = true;
		}
	}
	return;
}	

template<class T> inline void loadParamEnumText(bool& recognized, ifstream& file, string& id, const char* tag, T& target, string* compareList, T compareCount)
{
	if (!recognized)
	{
		if (id == tag)
		{
			int i = 0;
			string temp;
			file >> temp;
			while(i<compareCount)
			{
				if (temp == compareList[i])
				{
					target = (T)i;
					break;
				}
				i++;
			}
			if (i == compareCount)
			{
				cerr << "ERROR: No valid value for " << tag << ": " << temp << "\n";
				exit(-1);
			}

			cout << "Set " << tag << " : " << target << " (" << temp << ")\n";
			recognized = true;
		}
	}
	return;
}	


	
bool Parameters::readFromFile(const char* pf, bool initialRun)
{
	if ((pf == NULL) || (*pf == 0))
		return false;

 	ifstream parFile;

	parFile.open(pf);
 	if (!parFile.good())
	{
		cerr << "Could not open parameter file [" << pf << "].\n";
		return false;
	}

	cout << "Reading parameter file [" << pf << "]\n";

	string id;
  	bool recognized;

    while (parFile.good() )
    {
		parFile >> id;
		if (!parFile.good()) break;
      	if (id[0] == '#')  // if not a comment....
      		getline(parFile,id);
      	else
      	{
	      	recognized = false;

			// This prevents unwanted new seeding upon loading new parameters. 
			// preSeededSeedDensity has to be restated explicitly in the new parameter file if new seeds are desired.
			preSeededSeedDensity = 0.; 
			

			// the following parameters are only read upon program initialization
			// they cannot be changed by on-the-fly parameter updates
			if (initialRun)
			{
				loadParam(recognized, parFile, id, "discreteAngleNumber", discreteAngleNumber);
				loadParamRandom(recognized, parFile, id, "random_seed", seed);	// special

				if (!recognized)
				{
					loadParamEnumText<GeometryType>(recognized, parFile, id, "geometry", geometry, GeometryTypeText, g_COUNT_LAST);
					if (recognized)
					{
						// read the rest of the parameters
						parFile >> geomParam1;
						cout << "\tGeometry parameter 1 : " << geomParam1 << "\n";
						parFile >> geomParam2;
						cout << "\tGeometry parameter 2 : " << geomParam2 << "\n";
						parFile >> geomParam3;
						cout << "\tGeometry parameter 3 : " << geomParam3 << "\n";
					}

				}
			}

			loadParam(recognized, parFile, id, "outputDir", outputDir);
			loadParam(recognized, parFile, id, "createSubdir", createSubdir);

			loadParam(recognized, parFile, id, "newParameterReadInterval", newParameterReadInterval);
			loadParam(recognized, parFile, id, "newParameterFile", newParameterFile);

			loadParam(recognized, parFile, id, "zippering", zipperingEnabled);
			loadParam(recognized, parFile, id, "ind_cat", catastrophesEnabled);

			loadParam(recognized, parFile, id, "proportionalCatastrophes", proportionalCatastrophes);
			loadParam(recognized, parFile, id, "induced_cat_fraction", inducedCatastropheFraction);
			loadParam(recognized, parFile, id, "cat_start_angle", catStartAngle);
			loadParam(recognized, parFile, id, "zipFraction", zipFraction);
			loadParam(recognized, parFile, id, "magic_angle", magicAngle);
			loadParam(recognized, parFile, id, "c0Value", c0Value);
			loadParam(recognized, parFile, id, "z0Value", z0Value);


			loadParam(recognized, parFile, id, "vPlus", vPlus);
			loadParam(recognized, parFile, id, "vMin", vMin);
			loadParam(recognized, parFile, id, "vTM", vTM);
			
			loadParam(recognized, parFile, id, "numberOfRays", numberOfRays);
			loadParam(recognized, parFile, id, "maxMetaTrajectoryLength", maxMetaTrajectoryLength);
			loadParam(recognized, parFile, id, "diffusionCoefficient", diffusionCoefficient);
			loadParam(recognized, parFile, id, "rvRejectionUnbound", rvRejectionUnbound);
			loadParam(recognized, parFile, id, "rvRejectionBound", rvRejectionBound);
			loadParam(recognized, parFile, id, "unboundNucleationRate", unboundNucleationRate);
			
#ifdef BAND_CAT
      if (! recognized) {
        loadParam(recognized, parFile, id, "kCat", kCat[0]);
        if (recognized) {
          parFile >> kCat[1];
          cout << "\tSecond kCat value : " << kCat[1] << "\n";
        }
      }
      if (! recognized) {
        loadParam(recognized, parFile, id, "bandGapWidth", bandGapWidth[0]);
        if (recognized) {
          parFile >> bandGapWidth[1];
          cout << "\tGap width : " << bandGapWidth[1] << "\n";
        }
      }
      loadParamEnumText<OrientationType>(recognized, parFile, id, "kCatOrient", kCatOrient, OrientationTypeText, o_COUNT_LAST);
			loadParam(recognized, parFile, id, "nSpirals", nSpirals);
      loadParamEnumText<PatternType>(recognized, parFile, id, "catPattern", catPattern, PatternTypeText, p_COUNT_LAST);
			loadParam(recognized, parFile, id, "calcNucBiasAngle", calcNucBiasAngle);
			loadParam(recognized, parFile, id, "redistributeNucleationOverBands", redistributeNucleationOverBands);
			loadParam(recognized, parFile, id, "reducedGapNucleation", reducedGapNucleation);
			loadParam(recognized, parFile, id, "gapNucleationAcceptFraction", gapNucleationAcceptFraction);
#else
			loadParam(recognized, parFile, id, "kCat", kCat);
#endif
			loadParam(recognized, parFile, id, "kRes", kRes);
			loadParam(recognized, parFile, id, "kResExtraMax", kResExtraMax);
			loadParam(recognized, parFile, id, "kResExtraAngle", kResExtraAngle);
			loadParamEnumText<extraRescueType>(recognized, parFile, id, "extraRescueFunction", extraRescueFunction, extraRescueTypeText, r_COUNT_LAST);
			loadParam(recognized, parFile, id, "kSev", kSev);
			loadParam(recognized, parFile, id, "kCross", kCross);
			loadParam(recognized, parFile, id, "poolDensity", poolDensity);
			loadParam(recognized, parFile, id, "treadmilling", treadmillingEnabled);
			loadParam(recognized, parFile, id, "severing", severingEnabled);
			loadParam(recognized, parFile, id, "crossSevering", crossSeveringEnabled);
			loadParam(recognized, parFile, id, "crossSeveringTop", crossSeveringTop);
			loadParam(recognized, parFile, id, "crossSeveringTopFraction", crossSeveringTopFraction);
			loadParam(recognized, parFile, id, "crossSeveringStartAngle", crossSeveringStartAngle);
			loadParam(recognized, parFile, id, "catastropheMultiplier", catastropheMultiplier);
			loadParam(recognized, parFile, id, "forbiddenZones", forbiddenZones);
			loadParam(recognized, parFile, id, "restrictedPool", restrictedPool);

			loadParam(recognized, parFile, id, "edgeCatastropheEnabled", edgeCatastropheEnabled);
			loadParam(recognized, parFile, id, "edgeCatastropheSmooth", edgeCatastropheSmooth);
			loadParam(recognized, parFile, id, "pCatRegularEdge", pCatRegularEdge);
			loadParam(recognized, parFile, id, "pCatSpecialEdge", pCatSpecialEdge);
#ifdef VAR_CAT
			loadParam(recognized, parFile, id, "kCatAlpha", kCatAlpha);
			loadParam(recognized, parFile, id, "kCatBands", kCatBands);
			loadParamEnumText<OrientationType>(recognized, parFile, id, "kCatOrient", kCatOrient, OrientationTypeText, o_COUNT_LAST);
#endif
			loadParam(recognized, parFile, id, "kNuc", kNuc);

			loadParam(recognized, parFile, id, "stopTime", stopTime);
			loadParam(recognized, parFile, id, "measurementInterval", measurementInterval);
			loadParam(recognized, parFile, id, "densityLimit", densityLimit);
			loadParam(recognized, parFile, id, "wallClockLimit", wallClockLimit);
			loadParam(recognized, parFile, id, "memoryLimit", memoryLimit);
			loadParam(recognized, parFile, id, "movieEnabled", movieEnabled);
			loadParam(recognized, parFile, id, "movieFrameInterval", movieFrameInterval);
			loadParam(recognized, parFile, id, "angleHistogramBins", angleHistogramBins);
                        loadParam(recognized, parFile, id, "angleHistogramOptical", angleHistogramOptical);

			loadParam(recognized, parFile, id, "hiresLengthHistogramBins", hiresLengthHistogramBins);
			loadParam(recognized, parFile, id, "loresLengthHistogramBins", loresLengthHistogramBins);
			loadParam(recognized, parFile, id, "loresAngleHistogramBins", loresAngleHistogramBins);
			loadParam(recognized, parFile, id, "loresLifetimeHistogramBins", loresLifetimeHistogramBins);
			loadParam(recognized, parFile, id, "hiresLifetimeHistogramBins", hiresLifetimeHistogramBins);
			loadParam(recognized, parFile, id, "histogramAverageSamples", histogramAverageSamples);
			loadParamEnumText<SpatialHistogramType>(recognized, parFile, id, "spatialHistogramType", spatialHistogramType, SpatialHistogramTypeText, sh_COUNT_LAST);
			loadParam(recognized, parFile, id, "spatialHistogramBinsX", spatialHistogramBinsX);
			loadParam(recognized, parFile, id, "spatialHistogramBinsY", spatialHistogramBinsY);
			loadParam(recognized, parFile, id, "spatialHistogramCountCaps", spatialHistogramCountCaps);
			loadParam(recognized, parFile, id, "spatialHistogramWriteRaw", spatialHistogramWriteRaw);
			loadParam(recognized, parFile, id, "regionalOutputQuantities", regionalOutputQuantities);

			loadParamEnumText<NucleationType>(recognized, parFile, id, "nucleationType", nucleationType, NucleationTypeText, nuc_COUNT_LAST);
			loadParamEnumText<NucleationBiasType>(recognized, parFile, id, "nucleationBiasType", nucleationBiasType, NucleationBiasTypeText, nbias_COUNT_LAST);
			loadParam(recognized, parFile, id, "nucleationAlpha", nucleationAlpha);
			loadParam(recognized, parFile, id, "nucleationBiasAngle", nucleationBiasAngle);


			loadParamEnumText<InteractionType>(recognized, parFile, id, "interactionType", interactionType, InteractionTypeText, int_COUNT_LAST);
			loadParamEnumText<BundleType>(recognized, parFile, id, "bundleType", bundleType, BundleTypeText, bdl_COUNT_LAST);


			loadParam(recognized, parFile, id, "nucleationHalfIsotropicDensity", nucleationHalfIsotropicDensity);
                        loadParam(recognized, parFile, id, "ellipseReducedFreeRate", ellipseReducedFreeRate);
                        loadParam(recognized, parFile, id, "ellipseReducedFreeRateAcceptFraction", ellipseReducedFreeRateAcceptFraction);
			loadParam(recognized, parFile, id, "ellipseEpsilon", ellipseEpsilon);
			loadParam(recognized, parFile, id, "ellipseEpsilonAlongMT", ellipseEpsilonAlongMT);
			loadParam(recognized, parFile, id, "ellipseForwardAlongMT", ellipseForwardAlongMT );
			loadParam(recognized, parFile, id, "ellipseLeftFraction", ellipseLeftFraction);
			loadParam(recognized, parFile, id, "ellipseRightFraction", ellipseRightFraction);
			loadParam(recognized, parFile, id, "ellipseBackwardFraction", ellipseBackwardFraction);
			loadParam(recognized, parFile, id, "ellipseSidewaysAngle", ellipseSidewaysAngle);
			loadParam(recognized, parFile, id, "ellipseHalfIsotropicDensity", nucleationHalfIsotropicDensity);
			loadParam(recognized, parFile, id, "useRegionalNucleationSaturation", useRegionalNucleationSaturation);
			loadParam(recognized, parFile, id, "f_unbound", f_unbound);
			loadParam(recognized, parFile, id, "useDoubleNucleationSaturation", useDoubleNucleationSaturation);
			loadParam(recognized, parFile, id, "rn_targetratio", rn_targetratio);
			loadParam(recognized, parFile, id, "occupancytimeNC", occupancytimeNC);
			loadParam(recognized, parFile, id, "gridAspectratio", gridAspectratio);
			loadParam(recognized, parFile, id, "useAbsoluteBidirectionalNucBias", useAbsoluteBidirectionalNucBias);
			loadParam(recognized, parFile, id, "absoluteNucBiasAngle", absoluteNucBiasAngle);
			loadParam(recognized, parFile, id, "nucBiasVariance", nucBiasVariance);
			loadParam(recognized, parFile, id, "useMTdeflection", useMTdeflection);
			loadParam(recognized, parFile, id, "deflectionStepsize", deflectionStepsize);
			loadParam(recognized, parFile, id, "deflectionMaxAngle", deflectionMaxAngle);
			loadParam(recognized, parFile, id, "deflectionMinAngle", deflectionMinAngle);
			loadParam(recognized, parFile, id, "deflectionBundleCompensation", deflectionBundleCompensation);
			loadParam(recognized, parFile, id, "deflectionBundleProtection", deflectionBundleProtection);
			loadParam(recognized, parFile, id, "biasDeflectionRight", biasDeflectionRight);
			loadParam(recognized, parFile, id, "useBundleTracking", useBundleTracking);
			loadParam(recognized, parFile, id, "bundleTrackMaxAngle", bundleTrackMaxAngle);
			loadParam(recognized, parFile, id, "outputNucPos", outputNucPos);

			loadParamEnumText<NucleationType>(recognized, parFile, id, "preSeededType", preSeededType, NucleationTypeText, nuc_COUNT_LAST);
			loadParam(recognized, parFile, id, "preSeededAlpha", preSeededAlpha);
			loadParam(recognized, parFile, id, "preSeededRate", preSeededRate);
			loadParam(recognized, parFile, id, "preSeededSeedDensity", preSeededSeedDensity); // This needs to be the very last entry in the parfile or things will go wrong!

	
			
			if (!recognized)
			{
				cerr << "Parameter not recognized: " << id << "\n";	 
			}
			getline(parFile, id);	// reads the remainder of the line
      	}
    }
    cout << "Finished loading parameters.\n";
	return true;
}

void Parameters::verifyParameters()                  //end in line 1275
{
	int i;
	
	vTM *= treadmillingEnabled;
	kSev *= severingEnabled;
	kCross *= crossSeveringEnabled;

	if (vTM >= vPlus)
	{
		cerr << "ERROR: treadmilling velocity larger than growth velocity. Microtubule growth is not possible. Exiting.\n";
		exit(-1);
	}

	if (vMin >= vPlus)
	{
		cerr << "ERROR: Shrinking microtubules overtake growing microtubules. This will lead to queueing problems. Exiting.\n";
		exit(-1);
	}

#ifdef BAND_CAT
	if ((vPlus <= 0) || (vTM <0) || (kCat[0]*kCat[1] < 0) || (kRes < 0) || (kNuc < 0) || (kSev < 0) || (kCross < 0) || (poolDensity < 0))
#else 
	if ((vPlus <= 0) || (vTM <0) || (kCat < 0) || (kRes < 0) || (kNuc < 0) || (kSev < 0) || (kCross < 0) || (poolDensity < 0))
#endif
	{
		cerr << "Only vMin is allowed to be negative. Exiting to avoid nonsensical results.\n";
		exit(-1);
	}
  if (angleHistogramOptical) {
    if (angleHistogramBins < 2) {
      cerr << "angleHistogramOptical is pointless with <2 bins (" << angleHistogramBins << ") Not calculating it.\n" ;
      angleHistogramOptical = 0;
    }
  }

	if (restrictedPool != 0)
	{	// if there is a finite tubulin pool, check the validity of poolDensity
		if ((poolDensity <= 0) || (poolDensity > VERY_LARGE))
		{
			cerr << "ERROR [finite tubulin pool]: maximum density is out of range (value="<< poolDensity << "\n";
			exit(-1);
		}
	}

	if (movieEnabled)
	{
		if (movieFrameInterval < ZERO_CUTOFF)
		{
			cerr << "ERROR: movie recording cannot be enabled without specifying 'movieFrameInterval'\n";
			exit(-1);
		}
	}

	switch (nucleationType)
	{
    case nuc_biased: 
      discreteAngleNumber = 0;
      if (abs(nucleationAlpha) < ZERO_CUTOFF) 
      {
        cout << "nucleationAlpha should not be zero for nucleation type 2\n Changing to nucleationType = 0 (uniform)\n";
        nucleationType = nuc_isotropic;
        nucleationAlpha = 0.0;
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// checks for nuc_aster
    case nuc_ellipse: // further checks for nuc_ellipse below!
      if (nucleationBiasType != nbias_cross) { // for maintenance: this is ellipsePolar and ellipseApolar. 
        if (nucleationBiasAngle < 0 || nucleationBiasAngle > 2*PI )
          cout << "Adjusting nucleationBiasAngle from " << nucleationBiasAngle;
        while (nucleationBiasAngle > 2*PI )
          nucleationBiasAngle -= 2*PI;
        while (nucleationBiasAngle < 0 )
          nucleationBiasAngle += 2*PI;
        cout << " to " << nucleationBiasAngle << "\n";
        if (nucleationAlpha < 0 || nucleationAlpha >= 1 ) {
          cerr << "ERROR: nucleationAlpha must be in [0,1) for nucleationBiasType ellipsePolar and ellipseApolar, but is: " <<  nucleationAlpha << " Exiting.\n";
          exit(-1);
        }
      }
      if (nucleationType == nuc_biased) {
        break; }
		//case nuc_ellipse: // moved up to have only one entry point; 
			if (nucleationHalfIsotropicDensity < ZERO_CUTOFF)
			{
				cout << "Ellipse nucleation: some background of isotropic nucleation is required (nucleationHalfIsotropicDensity > 0)!\n";
				exit(-33);
			}
			if (abs(ellipseLeftFraction - ellipseRightFraction) > ZERO_CUTOFF)
				cout << "WARNING! nucleation type ellipse with asymmetric parameters (LeftFraction != RightFraction)\r\n";
			if (ellipseLeftFraction + ellipseRightFraction + ellipseBackwardFraction > 1.)
			{
				cout << "Ellipse nucleation: ellipseLeftFraction + ellipseRightFraction + ellipseBackwardFraction > 1.\r\n";
				exit(-12);
			}
			discreteAngleNumber = 0;
			if (nucleationAlpha < 0 ) 
			{
				cout << "Ellipse nucleation: nucleationAlpha should be positive. Setting to 0.\n";
				nucleationAlpha = 0.0;
			}
      if (ellipseEpsilonAlongMT == -1)
      {
        ellipseEpsilonAlongMT = ellipseEpsilon;
      }
      if (ellipseEpsilonAlongMT == 1)
      {
        ellipseForwardAlongMT = 1;
      }
			if (!ellipseReducedFreeRate)
            	{ ellipseReducedFreeRateAcceptFraction = 1.;}
            else 
			{
				if (ellipseReducedFreeRateAcceptFraction < ZERO_CUTOFF)
				{
					cout << "WARNING! No initiation of system possible without some free nucleations.\n";
                }
                else if (ellipseReducedFreeRateAcceptFraction > 1. - ZERO_CUTOFF)
                {
					cout << "No actual rejection of free nucleation events. Switching it off. (It is impossible to increase the rate of free nucleations beyond that of bound nucleations.)\n";
                    ellipseReducedFreeRate = 0;
                    ellipseReducedFreeRateAcceptFraction = 1.;
                }
			}
			break;	

    case nuc_discreteAngles: 
      if(discreteAngleNumber < 1)
      {
        cout << "Number of discrete nucleation angles cannot be zero or negative. Exiting.\n";
        exit(-1);
      }
      else
			{
				nucleationAngles = new double[discreteAngleNumber];
				for (i=0; i<discreteAngleNumber; i++)
				{
					nucleationAngles[i] = PI*i/discreteAngleNumber;// - 0.5*PI; // all in [-1/2pi, 1/2pi] -- assuming no effect of direction.
				}
				cout << "Calculated nucleationAngles: ";
				for (i=0; i<discreteAngleNumber; i++)
				{
					cout << nucleationAngles[i] << '\t' ;
				}
				cout << "\n";
			}

			if (loresAngleHistogramBins > discreteAngleNumber)
			{
				loresAngleHistogramBins = discreteAngleNumber;
				cout << "Setting loresAngleHistogramBins to " << discreteAngleNumber << " [number of discrete nucleation angles]\n";
			}
			if (angleHistogramBins > discreteAngleNumber)
			{
				angleHistogramBins = discreteAngleNumber;
				cout << "Setting angleHistogramBins to " << discreteAngleNumber << " [number of discrete nucleation angles]\n";
			}



			break;

		case nuc_isotropic: 
			discreteAngleNumber = 0;
			break;
		case nuc_aster:
			break;
		default:
			cout << "Unknown nucleationType. Aborting.\n";
			exit(-1);
			break;	
	}

	if (nucleationType != nuc_discreteAngles)
	{
		discreteAngleNumber = 0;
	}
	if (nucleationType != nuc_biased && nucleationType != nuc_ellipse )
	{
		if (abs(nucleationAlpha) > ZERO_CUTOFF)
		{
			cout << "nucleationAlpha not used\n Set to 0.\n";
			nucleationAlpha = 0.0;
		}
	}

	if (preSeededSeedDensity < -ZERO_CUTOFF)
	{
		cerr << "preSeededSeedDensity cannot be negative. Set to 0.\n";
		preSeededSeedDensity = 0;
	}
	if (preSeededRate < -ZERO_CUTOFF)
	{
		cerr << "preSeededRate cannot be negative. Set to 0.\n";
		preSeededRate = 0;
	}

	if (system->seedsLeft && preSeededSeedDensity > ZERO_CUTOFF)
	{
		cout << "NOTE: New nucleation seeds will replace old seeds. " << system->seedsLeft << " seeds destroyed.\n";
		system->seedsLeft=0;
	}

	if ((preSeededSeedDensity > ZERO_CUTOFF) && (preSeededRate < ZERO_CUTOFF))
	{
		cerr << "WARNING: New nucleation seeds will not be triggered because preSeededRate = 0.\n";
	}
	//if (preSeededSeedDensity > ZERO_CUTOFF) 
	//{
	//	if (preSeededRate > ZERO_CUTOFF)
	//	{
	//		double calcArea; // required because the geometry is created after verifying the initial parameter file. 
	//		switch (geometry) 
	//		{
	//			case g_periodic:
	//			case g_grid:
	//			case g_wormhole:
	//				calcArea = geomParam1 * geomParam2;
	//				break;
	//			case g_cylinder:
	//			case g_gridcylinder:
	//				if (forbiddenZones)
	//					calcArea = PI*(geomParam1*(2.*geomParam2));
	//				else
	//					calcArea = PI*(geomParam1*(geomParam1 + 2.*geomParam2));
	//				break;
	//			case g_box:
	//				if (forbiddenZones) 
	//				{
	//					// forbiddenZones: x,z plane not participating.
	//					calcArea = 2.*geomParam1*geomParam2 + 2.*geomParam2*geomParam3;
	//				}
	//				else
	//					calcArea = 2.*geomParam1 *( geomParam2 + geomParam3) + 2.*geomParam2*geomParam3;
	//				break;
	//			case g_pancake:
	//				calcArea = PI*geomParam1*geomParam1;
	//				break;
	//			default:
	//				cerr << "ERROR: no valid geometry handler present. aborting.\n";
	//				exit(-1);
	//		} 
	//		system->seedsLeft += static_cast<int> (0.499 + preSeededSeedDensity * calcArea) ;
	//		cout << "preSeeding "<< system->seedsLeft << " special nucleation events.\n";
	//	}
	//	else 
	//	{
	//		cout << "Trying to put seeds for nucleation that will never fire. Exit.\n";
	//		exit(-43);
 //   
	//	}
 // 
	//}

	if (preSeededSeedDensity > ZERO_CUTOFF) 
	{
		switch(preSeededType)
		{
			case nuc_isotropic:
				if (preSeededAlpha > ZERO_CUTOFF) 
				{
					cout << "isotropic pre-seeded nucleation implies preSeededAlpha = 0. Set to 0.\n";
					preSeededAlpha = 0.;
				}
				break;
			case nuc_biased:
				if (abs(preSeededAlpha) < ZERO_CUTOFF) 
				{
					cout << "biased pre-seeded nucleation with |alpha| = 0 is actually isotropic. Set to isotropic.\n";
					preSeededType = nuc_isotropic;
				}
				break; 
			case nuc_aster:
			case nuc_ellipse:
			case nuc_discreteAngles:
			default:
				cerr << "Invalid type for preSeeded nucleations: "<< preSeededType << ". Exit.\n";
				exit(-44);
		}
	}

	
	switch(interactionType)
  {
    case int_catFirst: // nothing special: no corrections
      break;
    case int_minimalFourier:
      if (!catastrophesEnabled)
      {
        cout << "minimalFourier interaction mode requires catastrophic collisions. Catastrophes enabled.\n";
        catastrophesEnabled = true;
      }
      if ((inducedCatastropheFraction < ZERO_CUTOFF) || (inducedCatastropheFraction > 1.0))
      {
        cout << "catastropheFraction must be positive and not larger than 1. Set to 1.\n";
        inducedCatastropheFraction = 1.0;
      }
      if (((c0Value - 0.75) < -ZERO_CUTOFF) || ((c0Value - 1.125) > ZERO_CUTOFF))
      {
        c0Value = 1.0;
        cout << "c0Value must be in the range [3/4, 9/8]. Set to 1.\n";
		}
		if (zipperingEnabled && (z0Value  < ZERO_CUTOFF))
		{
			zipperingEnabled = false;
			z0Value = 0;
			cout << "z0Value too low. Zippering disabled.\n";
		}

		break;
	case int_zipFirst:
		if (zipperingEnabled){
			if (zipFraction < 0.)
				zipFraction = 0.;
			if (zipFraction > 1.)
				zipFraction = 1.;
		}
		break;
	}

	if (edgeCatastropheEnabled)
	{
		if (pCatRegularEdge < -ZERO_CUTOFF)
			pCatRegularEdge = 0.0;
		if (pCatSpecialEdge < -ZERO_CUTOFF)
			pCatSpecialEdge = 0.0;
		if (pCatRegularEdge > 1.0)
			pCatRegularEdge = 1.0;
		if (pCatSpecialEdge > 1.0)
			pCatSpecialEdge = 1.0;

	}

	if (densityLimit <= 0)
	{
		cerr << "densityLimit should be a positive number. Ignoring limit.\n";
		densityLimit = VERY_LARGE;
	}

#ifdef VAR_CAT
    catMin = kCat * (1. - kCatAlpha);
    catMax = kCat * (1. + kCatAlpha);
  if (kCatBands > 0) {
    if ( geometry == g_sphere || geometry == g_pancake ) {
      cerr << "Variable kCat not implemented for geometry type " << GeometryTypeText[geometry] << ".\n";
      exit(-235);
    } 
    //catMin = kCat * (1. - kCatAlpha);
    Gmax = calcG(catMin,kRes);
    if (Gmax >= 0 && restrictedPool == 0) {
      cerr << "Entering unbound growth regime! (kCatMin: " << catMin << ", G0max: " << Gmax << ").\n";
      exit(-124);
    }
    cout <<  "kCatMin: " << catMin << ", G0max: " << Gmax << ".\n";
    //catMax = kCat * (1. + kCatAlpha);
    Gmin = calcG(catMax,kRes);
    cout <<  "kCatMax: " << catMax << ", G0min: " << Gmin << ".\n";
    switch (kCatOrient) {
      case o_x: 
        invKCatPeriod = 2.*PI*kCatBands / geomParam1; // band around x-axis. 
        break;
      case o_y:
        invKCatPeriod = 2.*PI*kCatBands / geomParam2; // band around y-axis. 
        break;
      case o_z:
        if (geometry == g_box) {
          invKCatPeriod = 2.*PI*kCatBands / geomParam3; // band around z-axis. 
        }
        else {
          cerr << "Invalid orientation for geometry type "<<  GeometryTypeText[geometry] << ".\n";
          exit (-125);
        }
        break;
    }
  }
#endif

#ifdef BAND_CAT
    catMin = min(kCat[0], kCat[1]);
    catMax = max(kCat[0], kCat[1]);
  if (bandGapWidth[0] > 0.) {
    if ( geometry == g_sphere || geometry == g_pancake ) {
      cerr << "Variable kCat not implemented for geometry type " << GeometryTypeText[geometry] << ".\n";
      exit(-235);
    } 
    //catMin = min(kCat[0], kCat[1]);
    Gmax = calcG(catMin,kRes);
    if (Gmax >= 0 && restrictedPool == 0) {
      cerr << "Entering unbound growth regime! (kCatMin: " << catMin << ", G0max: " << Gmax << ").\n";
      exit(-124);
    }
    cout <<  "kCatMin: " << catMin << ", G0max: " << Gmax << ".\n";
    //catMax = max(kCat[0], kCat[1]);
    Gmin = calcG(catMax,kRes);
    cout <<  "kCatMax: " << catMax << ", G0min: " << Gmin << ".\n";

    if (nSpirals > 0) {
      switch (kCatOrient) {
        case o_x: 
          switch (geometry){
            case g_grid: 
            case g_periodic: 
              wrapLength = geomParam2;
              //projectedPeriod = geomParam2*atan(asin((bandGapWidth[0]+bandGapWidth[1])/geomParam2)); // correct for angle: width of band and band-band distance are input; these are internal numbers.
              ////spiralPitch = nSpirals *(bandGapWidth[0]+bandGapWidth[1]) /(geomParam2); // band around x-axis. 
              //spiralPitch = nSpirals *(projectedPeriod) /(geomParam2); // band around x-axis. 
              break;
            case g_cylinder:
            case g_gridcylinder:
              wrapLength = 2.*PI*geomParam2;
              //wrapLength = geomParam1;  //TEMP
              //projectedPeriod = 2.*PI*geomParam2*atan(asin((bandGapWidth[0]+bandGapWidth[1])/(2.*PI*geomParam2)));
              ////spiralPitch = nSpirals *(bandGapWidth[0]+bandGapWidth[1]) /(2.*PI*geomParam2); // band around x-axis. 
              //spiralPitch = nSpirals *(projectedPeriod) /(2.*PI*geomParam2); // band around x-axis. 
              break;
            case g_box: 
              wrapLength = 2.*(geomParam2+geomParam3);
              //projectedPeriod = 2.*(geomParam2+geomParam3)*atan(asin((bandGapWidth[0]+bandGapWidth[1])/(2.*(geomParam2+geomParam3))));
              //spiralPitch = nSpirals *(projectedPeriod) /(2.*(geomParam3+geomParam2)); // band around x-axis. 
              break;
          }
          break;
        case o_y:
          switch (geometry){
            case g_grid: 
            case g_periodic: 
              wrapLength = geomParam1;
              break;
            case g_cylinder:
            case g_gridcylinder:
              cout << "Orientation y for spiral on cylinder is problematic.\n";
              wrapLength = geomParam1;
              //wrapLength = 2.*PI*geomParam2; //TEMP
              //projectedPeriod = geomParam1*atan(asin((bandGapWidth[0]+bandGapWidth[1])/geomParam1)); // correct for angle: width of band and band-band distance are input; these are internal numbers.
              //spiralPitch = nSpirals *(projectedPeriod) /(geomParam1); // band around y-axis. 
              break;
            case g_box: 
              wrapLength = 2.*(geomParam1+geomParam3);
              //projectedPeriod = 2.*(geomParam1+geomParam3)*atan(asin((bandGapWidth[0]+bandGapWidth[1])/(2.*(geomParam1+geomParam3))));
              //spiralPitch = nSpirals *(projectedPeriod) /(2.*(geomParam1+geomParam3)); // band around y-axis. 
              break;
          }
          break;
        case o_z:
          if (geometry == g_box) {
              wrapLength = 2.*(geomParam1+geomParam2);
            //projectedPeriod = 2.*(geomParam1+geomParam2)*atan(asin((bandGapWidth[0]+bandGapWidth[1])/(2.*(geomParam1+geomParam2))));
            //spiralPitch = nSpirals *(projectedPeriod) /(2.*(geomParam1+geomParam2)); // band around z-axis. 
          }
          else {
            cerr << "Invalid orientation for geometry type "<<  GeometryTypeText[geometry] << ".\n";
            exit (-125);
          }
          break;
      }

      projectedPeriod = wrapLength*atan(asin((bandGapWidth[0]+bandGapWidth[1])/wrapLength)); // correct for angle: width of band and band-band distance are input; these are internal numbers.
      spiralPitch = nSpirals *(projectedPeriod) /wrapLength; 
      projectedBand = bandGapWidth[0]/(bandGapWidth[0]+bandGapWidth[1])*projectedPeriod;
    }

    else {
      wrapLength = 1.; // irrelvant if nSpirals == 0
      spiralPitch = 0.;
      projectedPeriod = bandGapWidth[0]+bandGapWidth[1]; // no need to do projection in bands; use input values.
      projectedBand = bandGapWidth[0]; // no need to do projection in bands; use input values.
    }
    if (calcNucBiasAngle) {  // CHECK CALCULATION!
      nucleationBiasAngle = -atan(spiralPitch);
      if (kCatOrient == o_x) // Check what to do with o_z! ; anyway can't be uniform on a box!
        nucleationBiasAngle = 0.5*PI + nucleationBiasAngle ; 
      while (nucleationBiasAngle < 0 )  
        nucleationBiasAngle += 2.*PI; 
      while (nucleationBiasAngle > 2*PI )  
        nucleationBiasAngle -= 2.*PI; 
    }
    if (redistributeNucleationOverBands) {
      if (nucleationType == nuc_ellipse && ellipseForwardAlongMT ) { // Expect very messy behaviour if ellipseForwardAlongMT==1, because the same parent might trajectory induce multiple events in a neighbouring band, which will result in many overlapping trajectories.
        cerr << "redistributeNucleationOverBands does not work in combination with ellipseForwardAlongMT == 1. Exit.\n";
        exit(-126);
      }
      if ( geometry == g_box || geometry == g_pancake || geometry == g_sphere ) {
        cerr << "Invalid geometry for redistributeNucleationOverBands "<<  GeometryTypeText[geometry] << ".\n";
        exit (-125);
      }
      if (nSpirals == 0) 
        redistributeShift = bandGapWidth[0]+bandGapWidth[1] ;
      else 
        redistributeShift = projectedPeriod / nSpirals;
      switch (kCatOrient) {
        case o_x: 
          redistributeShiftNumber = geomParam1/redistributeShift;
          redistributeShiftMax = geomParam1/2.;
          break;
        case o_y: 
          redistributeShiftNumber = geomParam2/redistributeShift;
          redistributeShiftMax = geomParam2/2. ;
          break;
        case o_z: 
          cerr << "Invalid orientation in combination with redistributeNucleationOverBands "<<  GeometryTypeText[geometry] << ".\n";
          exit (-125);
          break;
      }
    }
  }

  // Test if max kRes obeys G<0
  if (kResExtraMax > 0) {
  //if (0) {
	    double resMax = kRes + kResExtraMax;
	    Gmax = calcG(catMin,resMax);
	    if (Gmax >= 0 && restrictedPool == 0) {
		    cerr << "Entering unbound growth regime! (kResMax: " << kRes+kResExtraMax << ", G0max: " << Gmax << ").\n";
		    exit(-124);
	    }
	    cout <<  "kResMax: " << kRes+kResExtraMax << ", G0max maximum: " << Gmax << ".\n";
    }
#else
  if (kResExtraMax > 0) {
    Gmax = calcG(kCat,kRes+kResExtraMax);
    if (Gmax >= 0 && restrictedPool == 0) {
      cerr << "Entering unbound growth regime! (kResMax: " << kRes+kResExtraMax << ", G0max: " << Gmax << ").\n";
      exit(-124);
    }
    cout <<  "kResMax: " << kRes+kResExtraMax << ", G0max: " << Gmax << ".\n";
  }
#endif

  if (spatialHistogramType != sh_none) {
    if (geometry != g_grid && geometry != g_gridcylinder ){
      cerr << "spatialMeasurement not implemented for geometry type " << GeometryTypeText[geometry] << ".\n";
      exit (-127);
    }
    if (spatialHistogramBinsX <= 0 && ( spatialHistogramType == sh_xy || spatialHistogramType == sh_x )) {
      cerr << "invalid number of bins for spatialHistogram in X-direction: " << spatialHistogramBinsX << ".\n";
      exit (-124);
    }
    if (spatialHistogramBinsY <= 0 && ( spatialHistogramType == sh_xy || spatialHistogramType == sh_y )) {
      cerr << "invalid number of bins for spatialHistogram in Y-direction: " << spatialHistogramBinsY << ".\n";
      exit (-1243);
    }
  }
  return;
}




bool Parameters::calcTheoryParameters(void)
{
	double maRad = magicAngle*PI/180.0;
	double csRad = catStartAngle*PI/180.0;

	switch (interactionType)
	{
		case int_zipFirst:
			if (zipperingEnabled)
				z0calc = (4./PI)*(1. - cos(maRad))*zipFraction;
			else
				z0calc = 0;

			if (catastrophesEnabled)
			{
				if (proportionalCatastrophes)
				{
					if (magicAngle > catStartAngle)
						c0calc = (4.0/PI)*2.0*inducedCatastropheFraction*(1.- sin(maRad)+(maRad-csRad)*cos(maRad))/(PI - 2.*csRad);
					else
						c0calc = (4.0/PI)*2.0*inducedCatastropheFraction*(1.0 - sin(csRad))/(PI - 2.0*csRad);
				}
				else
				{
					c0calc = (4.0/PI)*inducedCatastropheFraction*cos(max(maRad,csRad));
				}
			}
			else
				c0calc = 0;
			x0calc = (4.0/PI) - z0calc - c0calc;
			break;
		case int_catFirst:

			if (catastrophesEnabled)
			{
				if (proportionalCatastrophes)
				{
					double temp=max(0.,csRad);
					c0calc = (4.0/PI)*2.0*inducedCatastropheFraction*(1.0-sin(temp)+(temp - csRad)*cos(temp))/(PI-2.0*csRad);
				}
				else
				{
					c0calc = (4.0/PI)*inducedCatastropheFraction*cos(max(csRad,0.));
				}
			}

			if (zipperingEnabled)
			{
				z0calc = (4.0/PI)*(1.0 - cos(maRad));
				if (catastrophesEnabled && (csRad < maRad))
				{
					if (proportionalCatastrophes)
					{
						double leftBoundary = max(0.,csRad);
						z0calc -= (4.0/PI)*2.0*inducedCatastropheFraction*((leftBoundary - csRad)*cos(leftBoundary) \
							- (maRad - csRad)*cos(maRad) - sin(leftBoundary) + sin(maRad))/(PI- 2.0*csRad);
					}
					else
					{
						z0calc -= (4.0/PI)*inducedCatastropheFraction*(cos(csRad)-cos(maRad)) ;
					}
				}
			}
			else
				z0calc = 0;

			x0calc = (4.0/PI) - z0calc - c0calc;
			break;
		case int_minimalFourier:
			c0calc = c0Value * inducedCatastropheFraction;
			z0calc = z0Value * inducedCatastropheFraction;
			x0calc = 4./PI - c0calc - z0calc;

			break;
		default:
			cerr << "no function to calculate z0, c0, x0 etc. (interactionType = " << interactionType << ")\n";
			z0calc = 123;
			c0calc = 123;
			x0calc = 123;
			return false ;

	}

#ifdef BAND_CAT
	double littleG = kRes/(-vMin+vTM) - kCat[0]/(vPlus - vTM);
#else
	double littleG = kRes/(-vMin+vTM) - kCat/(vPlus - vTM);
#endif
	lengthFactor = pow((2.0*(-vMin+vTM)*(vPlus-vTM))/(kNuc*(-vMin+vPlus)), 1./3.);
	G_0 = littleG * lengthFactor ;
	adjustedG = G_0 * pow((vPlus - vTM)/(vPlus),1./3.);

	cout << "\n\n# Calculated theory parameters \r\n";
	cout << "c0\t" << c0calc << "\r\n";
	cout << "z0\t" << z0calc << "\r\n";
	cout << "x0\t" << x0calc << "\r\n";
#ifdef BAND_CAT
	cout << "l_free (band)\t" << -1./littleG << " micrometer\r\n";
	cout << "G_0 (band)\t" << G_0 << "\r\n";
	cout << "adjustedG (band)\t" << adjustedG << "\r\n";
#else
	cout << "l_free\t" << -1./littleG << " micrometer\r\n";
	cout << "G_0\t" << G_0 << "\r\n";
	cout << "adjustedG\t" << adjustedG << "\r\n";
#endif
	cout << "l0\t" << lengthFactor << "\r\n\n";

  if (useMTdeflection)
  {
    double pl = - deflectionStepsize / ( 2*log( (sin(deflectionMaxAngle)-sin(deflectionMinAngle)+deflectionMinAngle)/deflectionMaxAngle ) );
    cout << "Persistence length: " << pl << "um\r\n\n";
  }

	return true;
}

double Parameters::calcG(double rCat , double rRes) { // this is lazy implementation: only kCat and kRes (locally called rCat) serve as input parameters.
	double littleG = rRes/(-vMin+vTM) - rCat/(vPlus - vTM);
	double lengthFactorAdjusted = pow((2.0*(-vMin+vTM)*pow((vPlus-vTM),2.))/(kNuc*vPlus*(-vMin+vPlus)), 1./3.);
	cout << "kCat: " << rCat << "; kRes: " << rRes <<'\n';
	return  littleG * lengthFactorAdjusted ;
}
