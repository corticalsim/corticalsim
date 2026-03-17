#include "corticalSimReal.h"

string DirectionTypeText[] = {"forward","backward"};
string GeometryTypeText[] = {"cell"};
string InteractionTypeText[] = {"zipFirst","catFirst","minimalFourier"};
string NucleationTypeText[] = {"isotropic","biased","ellipse" };
string BundleTypeText[] = {"simple","sticky","noZip","multiCollision","Ncollision" };
Parameters::Parameters(System* s) :
    system(s),
    vPlus(0.08),
    vMin(-0.16),
    vTM(0.01),
    kSev(0.001),
    geometryScaleFactor(15.0),
    kCross(0.001),
    kCat(0.005),
    kRes(0.007),
    poolDensity(5.0),
    treadmillingEnabled(1),
    severingEnabled(0),
    crossSeveringEnabled(0),
    edgeCatastropheEnabled(0),
    edgeCatastropheSmooth(0),
    pCatSpecialEdgeMax(0.0),
    pCatRegularEdgeMax(0.0),
    edgNumber(0),
    faceNumber(0),
    kNuc(0.001),
    nucleationType(nuc_isotropic),
    interactionType(int_zipFirst),
    zipperingEnabled(1),
    catastrophesEnabled(1),
    PPB(0),
    PPBformingTime(VERY_LARGE),
    proportionalCatastrophes(0),
    inducedCatastropheFraction(0.5),
    catStartAngle(0),
    zipFraction(1.0),
    PPBkNucFraction(0.0),
    magicAngle(40),
    crossSeveringTop(1),
    crossSeveringStartAngle(0),
    restrictedPool(0),
    bundleType(bdl_simple),
    geomParam("0"),
    seed(time(0)),
    stopTime(36000),
    measurementInterval(600),
    wallClockLimit(1000),
    memoryLimit(10000000),
    outputDir("."),
    inputDir("."),
    createSubdir(1),
    newParameterReadInterval(VERY_LARGE),
    newParameterFile(""),
    movieEnabled(0),			
    movieFrameInterval(-1),
    geometry ("2D-plane_1_0_0"),
    showMesh(0),
    showOutput(0),
    c0calc(0),
    x0calc(0),
    z0calc(0)		
{
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
        return (true);
    }
    else
    return (false);
}

bool Parameters::writeToFile()
{
	ofstream & of = system->parameterFile;
	of << "# CorticalSim parameters\r\n";
	of << "########## Input-Output path ########################\r\n";
	of << "inputDir\t" << inputDir << "\r\n"; 	
	of << "outputDir\t" << outputDir << "\r\n"; 
	of << "newParameterReadInterval\t"<< newParameterReadInterval<< "\r\n"; 
	of << "createSubdir\t"<< createSubdir << "\r\n"; 
	of << "newParameterFile\t"<< newParameterFile << "\r\n"; 
	of << "########## Dynamic parameters ##########\r\n"; 
	of << "kSev\t"<< kSev	<< "\r\n"; 
	of << "kCross\t"<< kCross  << "\r\n"; 
	of << "vPlus\t"<< vPlus	<< "\r\n"; 
	of << "vMin\t"<< vMin	<< "\r\n"; 
	of << "vTM\t"<< vTM	<< "\r\n"; 
	of << "kCat\t"<<  kCat 	 << "\r\n"; 
	of << "kRes\t"<< kRes	<< "\r\n"; 
	of << "kNuc\t"<<  kNuc 	 << "\r\n"; 
	of << "magic_angle\t"<< magicAngle	<< "\r\n"; 
	of << "poolDensity\t"<< poolDensity	<< "\r\n"; 
	of << "nucleationType\t"<< nucleationType	<< "\r\n"; 
	of << "########## Zip-IndCat-Severing events ######\r\n"; 
	of << "treadmilling\t"<< treadmillingEnabled	<< "\r\n"; 
	of << "ind_cat\t"<< catastrophesEnabled	       << "\r\n"; 
	of << "cat_start_angle\t"<< catStartAngle	<< "\r\n"; 
	of << "induced_cat_fraction\t"<< inducedCatastropheFraction	<< "\r\n"; 
	of << "zippering\t"<< zipperingEnabled	<< "\r\n"; 
	of << "zipFraction\t"<< zipFraction    << "\r\n"; 
	of << "severing\t"<< severingEnabled	<< "\r\n"; 
	of << "crossSevering\t"<< crossSeveringEnabled	<< "\r\n"; 
	of << "crossSeveringTop\t"<< crossSeveringTop << "\r\n"; 
	of << "crossSeveringStartAngle\t"<< crossSeveringStartAngle << "\r\n"; 
	of << "########## Edge catastrophe specification ####\r\n";
	map<string, double>::iterator it;
	of << "########################\r\n";
	of << "edgNumber\t"<< edgNumber	<< "\r\n";
	for(it = edgCatMap.begin(); it!= edgCatMap.end(); it++)
	of << it->first <<"\t"<< it->second  << "\r\n";
	of << "########################\r\n";
	of << "faceNumber\t"<< faceNumber	<< "\r\n"; 
	for(it = faceCatMap.begin(); it!= faceCatMap.end(); it++)
	of << it->first <<"\t"<< it->second  << "\r\n"; 
	of << "########################\r\n";
	of << "edgeCatastropheSmooth\t"<< edgeCatastropheSmooth << "\r\n"; 
	of << "edgeCatastropheEnabled\t"<< edgeCatastropheEnabled << "\r\n";  
	of << "PPBkNucFraction\t"<< PPBkNucFraction << "\r\n"; 
	of << "proportionalCatastrophes\t"<< proportionalCatastrophes 	 << "\r\n"; 
	of << "########## specific conditions ###########\r\n"; 
	of << "restrictedPool\t"<< restrictedPool 	 << "\r\n";   
	of << "interactionType\t"<< interactionType << "\r\n"; 
	of << "bundleType\t"<< bundleType	<< "\r\n"; 
	of << "########## Random seed and run time #####\r\n"; 
	of << "memoryLimit\t"<< memoryLimit	<< "\r\n"; 
	of << "wallClockLimit\t"<< wallClockLimit	<< "\r\n"; 
	of << "stopTime\t"<< stopTime	<< "\r\n"; 
	of << "random_seed\t"<< seed	<< "\r\n"; 
	of << "########## Movie select ############\r\n"; 
	of << "movieEnabled\t"<< movieEnabled	<< "\r\n";  
	of << "movieFrameInterval\t"<< movieFrameInterval	<< "\r\n"; 
	of << "measurementInterval\t"<< measurementInterval	<< "\r\n"; 
	of << "########## Display output ##########\r\n"; 
	of << "showMesh\t"<<  showMesh 	<< "\r\n"; 
	of << "showOutput\t"<< showOutput 	 << "\r\n"; 
	of << "########## Geometry ################\r\n";  
	of << "geometryScaleFactor\t"<< geometryScaleFactor 	<< "\r\n"; 
	of << "geometry\t"<< geometry 	 << "\r\n"; 

    return true;
}

template<class T> inline void loadParam(bool& recognized, ifstream& file, string& id, const char* tag, T& target)
{
    if (!recognized)
    {
        if (id == tag)
        {
            file >> target;
            recognized = true;
        }
    }
    return;
}

template<class T> inline void loadParamMaps(bool& recognized, ifstream& file, string& id, const char* tag, T& target,map<string,T>& edgCatMap)
{
    if (!recognized)
    {
        if (id == tag)
        {
            file >> target;
	    string key = boost::lexical_cast<string>(tag);
	    edgCatMap[key] = target; 
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
                file >> temp;
                while (temp > boost::lexical_cast<string>(ULONG_MAX)) 
                {
                    cout << temp <<" Bigger than" <<boost::lexical_cast<string>(ULONG_MAX)<<endl;
                    temp = temp.substr(1);
                }
                target = boost::lexical_cast<T>(temp);
            }
            else
            {
                cout << "Using loadParamRandom for the wrong parameter: " << id << ", quitting.\n";
                exit(-1);
            }
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

    if(showOutput == 1)
    cout << "Reading parameter file [" << pf << "]\n";

    string id;
    bool recognized;
    while (parFile.good() )
    {
        parFile >> id;
        if (!parFile.good()) 
	break;

        if (id[0] == '#')
	getline(parFile,id);
        else
        {
            recognized = false;

            if (initialRun)
            loadParamRandom(recognized, parFile, id, "random_seed", seed); 
    
            loadParam(recognized, parFile, id, "outputDir", outputDir);
	    loadParam(recognized, parFile, id, "inputDir", inputDir);
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
            loadParam(recognized, parFile, id, "vPlus", vPlus);
            loadParam(recognized, parFile, id, "vMin", vMin);
            loadParam(recognized, parFile, id, "vTM", vTM);
            loadParam(recognized, parFile, id, "kCat", kCat);
            loadParam(recognized, parFile, id, "kRes", kRes);
            loadParam(recognized, parFile, id, "kSev", kSev);
            loadParam(recognized, parFile, id, "kCross", kCross);
            loadParam(recognized, parFile, id, "poolDensity", poolDensity);
            loadParam(recognized, parFile, id, "treadmilling", treadmillingEnabled);
            loadParam(recognized, parFile, id, "severing", severingEnabled);
            loadParam(recognized, parFile, id, "crossSevering", crossSeveringEnabled);
            loadParam(recognized, parFile, id, "crossSeveringTop", crossSeveringTop);
            loadParam(recognized, parFile, id, "crossSeveringStartAngle", crossSeveringStartAngle);
            loadParam(recognized, parFile, id, "restrictedPool", restrictedPool);
            loadParam(recognized, parFile, id, "edgeCatastropheEnabled", edgeCatastropheEnabled);
            loadParam(recognized, parFile, id, "edgeCatastropheSmooth", edgeCatastropheSmooth);
	    loadParam(recognized, parFile, id, "pCatSpecialEdgeMax", pCatSpecialEdgeMax);
	    loadParam(recognized, parFile, id, "pCatRegularEdgeMax", pCatRegularEdgeMax);
	    loadParam(recognized, parFile, id, "PPB", PPB);
	    loadParam(recognized, parFile, id, "PPBformingTime",PPBformingTime);
            loadParam(recognized, parFile, id, "kNuc", kNuc);
            loadParam(recognized, parFile, id, "stopTime", stopTime);
            loadParam(recognized, parFile, id, "measurementInterval", measurementInterval);
            loadParam(recognized, parFile, id, "wallClockLimit", wallClockLimit);
            loadParam(recognized, parFile, id, "memoryLimit", memoryLimit);
            loadParam(recognized, parFile, id, "movieEnabled", movieEnabled);
            loadParam(recognized, parFile, id, "movieFrameInterval", movieFrameInterval);
            loadParam(recognized, parFile, id, "geometryScaleFactor",geometryScaleFactor);
            loadParamEnumText<NucleationType>(recognized, parFile, id, "nucleationType", nucleationType, NucleationTypeText, nuc_COUNT_LAST);
            loadParamEnumText<InteractionType>(recognized, parFile, id, "interactionType", interactionType, InteractionTypeText, int_COUNT_LAST);
            loadParamEnumText<BundleType>(recognized, parFile, id, "bundleType", bundleType, BundleTypeText, bdl_COUNT_LAST);
	    loadParam(recognized, parFile, id, "showOutput", showOutput);
	    loadParam(recognized, parFile, id, "showMesh", showMesh);
	    loadParam(recognized, parFile, id, "geometry", geometry);
            loadParam(recognized, parFile, id, "edgNumber",edgNumber);

	    for(int edg = 1; edg <= edgNumber;edg++)
	    {
                double edgCatVal;
		string edgCatType1 = "edgCat_" + boost::lexical_cast<string>(edg);
		char *edgCatType = new char[edgCatType1.length() + 1];
		strcpy(edgCatType, edgCatType1.c_str());
	    	loadParamMaps(recognized, parFile,id,edgCatType,edgCatVal,edgCatMap);
		delete [] edgCatType;
	    }

            loadParam(recognized, parFile, id, "faceNumber",faceNumber);
	    for(int face = 1; face <= faceNumber;face++)
	    {
                double faceCatVal;
		string faceCatType1 = "faceCat_" + boost::lexical_cast<string>(face);
		char *faceCatType = new char[faceCatType1.length() + 1];
		strcpy(faceCatType, faceCatType1.c_str());
	    	loadParamMaps(recognized, parFile,id,faceCatType,faceCatVal,faceCatMap);
		delete [] faceCatType;
	    }
            if (!recognized)
            {
                cerr << "Parameter not recognized: " << id << "\n";
            }
            getline(parFile, id);	
        }
    }

    if(showOutput == 1)
    cout << "Finished loading parameters.\n";

    return true;
}

void Parameters::verifyParameters()
{
    int i;

    // adjust the parameter values in accordance with the swith (on/off) values
    vTM *= treadmillingEnabled;
    kSev *= severingEnabled;
    kCross *= crossSeveringEnabled;

    // treadmilling speed is not allowed to exceed the plus-end velocity
    if (vTM >= vPlus)
    {
        cerr << "ERROR: treadmilling velocity larger than growth velocity. Microtubule growth is not possible. Exiting.\n";
        exit(-1);
    }

    // minus-end velocity is not allowed to exceed the plus-end velocity
    if (vMin >= vPlus)
    {
        cerr << "ERROR: Shrinking microtubules overtake growing microtubules. This will lead to queueing problems. Exiting.\n";
        exit(-1);
    }

    // only minus-end velocity is allowed to be negative
    if ((vPlus <= 0) || (vTM <0) || (kCat < 0) || (kRes < 0) || (kNuc < 0) || (kSev < 0) || (kCross < 0) || (poolDensity < 0))
    {
        cerr << "Only vMin is allowed to be negative. Exiting to avoid nonsensical results.\n";
        exit(-1);
    }

    // if there is a finite tubulin pool, check the validity of poolDensity
    if (restrictedPool != 0)
    {
        if ((poolDensity <= 0) || (poolDensity > VERY_LARGE))
        {
            cerr << "ERROR [finite tubulin pool]: maximum density is out of range (0 : VERY_LARGE) :"<< poolDensity << "\n";
            exit(-1);
        }
    }

    // it is not possible to create movie without valid movie-interval
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
	// isotropic nucleation need no conditional check
    	case nuc_isotropic:
        break;

    	default:
        	cout << "Unknown nucleationType. Aborting.\n";
        	exit(-1);
        break;
    }

    switch(interactionType)
    {
        // catastrophe first, need no further conditional check
    	case int_catFirst: 
        break;

        // zip first needs further check
    	case int_zipFirst:

		// zippering fraction must be in the range: (0.0 - 1.0)
       		 if (zipperingEnabled)
        	{
            		if (zipFraction < 0.)
                	zipFraction = 0.;
            		if (zipFraction > 1.)
                	zipFraction = 1.;
        	}
        break;
    }

    return;
}

bool Parameters::calcTheoryParameters(void)
{
    // convert the angle values from degrees to radians
    double maRad = magicAngle*PI/180.0;
    double csRad = catStartAngle*PI/180.0;

    switch (interactionType)
    {
	// zippering first is requested   	
	case int_zipFirst:

		// fzippering is enabled
        	if (zipperingEnabled)
            	z0calc = (4./PI)*(1. - cos(maRad))*zipFraction;

		// zippering is disabled
        	else
            	z0calc = 0;

        	// catastrophe is enabled
        	if (catastrophesEnabled)
        	{
	    		// proportional catastrophe is enabled
            		if (proportionalCatastrophes)
            		{
				// magic angle > catastrophe start angle
                		if (magicAngle > catStartAngle)
                    		c0calc = (4.0/PI)*2.0*inducedCatastropheFraction*(1.- sin(maRad)+(maRad-csRad)*cos(maRad))/(PI - 2.*csRad);

				// magic angle <= catastrophe start angle
                		else
                    		c0calc = (4.0/PI)*2.0*inducedCatastropheFraction*(1.0 - sin(csRad))/(PI - 2.0*csRad);
            		}

	    		// proportional catastrophe is disabled 
            		else
			c0calc = (4.0/PI)*inducedCatastropheFraction*cos(max(maRad,csRad));           		
        	}

		// catastrophe is disabled
        	else
            	c0calc = 0;

		// calculate x0, using the values of c0 and/or z0
        	x0calc = (4.0/PI) - z0calc - c0calc;
        break;

    	// catastrophe first is requested  
    	case int_catFirst:

		// catastrophe is enabled
        	if (catastrophesEnabled)
        	{
			// proportional catastrophe is enabled
            		if (proportionalCatastrophes)
            		{
                		double temp=max(0.,csRad);
                		c0calc = (4.0/PI)*2.0*inducedCatastropheFraction*(1.0-sin(temp)+(temp - csRad)*cos(temp))/(PI-2.0*csRad);
            		}

			// proportional catastrophe is disabled 
            		else
                	c0calc = (4.0/PI)*inducedCatastropheFraction*cos(max(csRad,0.));
        	}

        	// zippering is enabled
		if (zipperingEnabled)
        	{
            		z0calc = (4.0/PI)*(1.0 - cos(maRad));

			// catastrophe is enabled and  (magic angle > catastrophe start angle)
            		if (catastrophesEnabled && (csRad < maRad))
            		{
				// proportional catastrophe is enabled
                		if (proportionalCatastrophes)
                		{
                    			double leftBoundary = max(0.,csRad);
                    			z0calc -= (4.0/PI)*2.0*inducedCatastropheFraction*((leftBoundary - csRad)*cos(leftBoundary)
                              		- (maRad - csRad)*cos(maRad) - sin(leftBoundary) + sin(maRad))/(PI- 2.0*csRad);
                		}

				// proportional catastrophe is disabled 
                		else                
                    		z0calc -= (4.0/PI)*inducedCatastropheFraction*(cos(csRad)-cos(maRad)) ;                
            		}
        	}

		// zippering is disabled
        	else
            	z0calc = 0;

        	// calculate x0, using the values of c0 and/or z0
        	x0calc = (4.0/PI) - z0calc - c0calc;
        break;

    default:
        cerr << "no function to calculate z0, c0, x0 etc. (interactionType = " << interactionType << ")\n";
        return false ;

    }

    // theoretial value of G
    double l0 = pow((2.0*(-vMin+vTM)*(vPlus-vTM)*(vPlus - vTM))/(kNuc*vPlus*(-vMin+vPlus)), 1./3.);
    double lavgInv = kRes/(-vMin+vTM) - kCat/(vPlus - vTM);
    double G_0 =  l0*lavgInv ;
 
    if(showOutput == 1)
    {
	    cout <<"==============================================="<<endl;
	    cout << " Calculated theory parameters"<<endl;
	    cout <<"==============================================="<<endl;
	    cout << "l0 = " << l0 <<endl;;
	    cout << "lavg = " << -1./lavgInv <<endl;
	    cout << "G_0 = -l0/lavg = " << G_0 <<endl;
	    cout <<"==============================================="<<endl;
   }

    return (true);
}
