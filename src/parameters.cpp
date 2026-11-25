#include <sstream>

#include "parameters.hpp"
#include "system.hpp"

string direction_type_text[] = { "forward", "backward" };
string geometry_type_text[] = { "cell" };
string interaction_type_text[] = { "zipFirst", "catFirst", "minimalFourier" };
string nucleation_type_text[] = { "isotropic", "biased", "ellipse" };
string bundle_type_text[] = { "simple", "sticky", "noZip", "multiCollision", "Ncollision" };

Parameters::Parameters(System* s):
    system(s),
    vPlus(0.08),
    vMin(-0.16),
    vTM(0.01),
    kSev(0.001),
    kCross(0.001),
    kCat(0.005),
    kRes(0.007),
    kNuc(0.001),
    PPBformingTime(VERY_LARGE),
    poolDensity(5.0),
    treadmillingEnabled(1),
    severingEnabled(0),
    crossSeveringEnabled(0),
    crossSeveringTop(1),
    crossSeveringStartAngle(0),
    PPB(0),
    restrictedPool(0),
    edgeCatastropheEnabled(0),
    edgeCatastropheSmooth(0),
    pCatSpecialEdgeMax(0.0),
    pCatRegularEdgeMax(0.0),
    edgNumber(0),
    faceNumber(0),
    nucleationType(NUC_ISOTROPIC),
    zipperingEnabled(1),
    catastrophesEnabled(1),
    proportionalCatastrophes(0),
    inducedCatastropheFraction(0.5),
    zipFraction(1.0),
    PPBkNucFraction(0.0),
    catStartAngle(0),
    magicAngle(40),
    interactionType(INT_ZIP_FIRST),
    bundleType(BDL_SIMPLE),
    geomParam("0"),
    seed(time(nullptr)),
    stopTime(36000),
    measurementInterval(600),
    wallClockLimit(1000),
    memoryLimit(10000000),
    inputDir("."),
    outputDir("."),
    createSubdir(1),
    newParameterReadInterval(VERY_LARGE),
    newParameterFile(""),
    movieEnabled(0),
    showOutput(0),
    showMesh(0),
    geometry("2D-plane_1_0_0"),
    movieFrameInterval(-1),
    geometryScaleFactor(15.0),
    c0calc(0),
    x0calc(0),
    z0calc(0)
{
    return;
}

void Parameters::initialize(const char* pf)
{
    if (!read_from_file(pf, true))
    {
        cout << "Failed to load parameters from file [" << pf << "]. Aborting\n";
        exit(-1);
    }
    verify_parameters();
    calc_theory_parameters();

    return;
}

bool Parameters::reinitialize(const char* pf)
{
    if (read_from_file(pf, false))
    {
        verify_parameters();
        calc_theory_parameters();
        return (true);
    }
    else
    {
        return (false);
    }
}

bool Parameters::write_to_file()
{
    ofstream& of = system->parameterFile;
    of << "# CorticalSim parameters\r\n";
    of << "########## Input-Output path ########################\r\n";
    of << "inputDir\t" << inputDir << "\r\n";
    of << "outputDir\t" << outputDir << "\r\n";
    of << "newParameterReadInterval\t" << newParameterReadInterval << "\r\n";
    of << "createSubdir\t" << createSubdir << "\r\n";
    of << "newParameterFile\t" << newParameterFile << "\r\n";
    of << "########## Dynamic parameters ##########\r\n";
    of << "kSev\t" << kSev << "\r\n";
    of << "kCross\t" << kCross << "\r\n";
    of << "vPlus\t" << vPlus << "\r\n";
    of << "vMin\t" << vMin << "\r\n";
    of << "vTM\t" << vTM << "\r\n";
    of << "kCat\t" << kCat << "\r\n";
    of << "kRes\t" << kRes << "\r\n";
    of << "kNuc\t" << kNuc << "\r\n";
    of << "magic_angle\t" << magicAngle << "\r\n";
    of << "poolDensity\t" << poolDensity << "\r\n";
    of << "nucleationType\t" << nucleationType << "\r\n";
    of << "########## Zip-IndCat-Severing events ######\r\n";
    of << "treadmilling\t" << treadmillingEnabled << "\r\n";
    of << "ind_cat\t" << catastrophesEnabled << "\r\n";
    of << "cat_start_angle\t" << catStartAngle << "\r\n";
    of << "induced_cat_fraction\t" << inducedCatastropheFraction << "\r\n";
    of << "zippering\t" << zipperingEnabled << "\r\n";
    of << "zipFraction\t" << zipFraction << "\r\n";
    of << "severing\t" << severingEnabled << "\r\n";
    of << "crossSevering\t" << crossSeveringEnabled << "\r\n";
    of << "crossSeveringTop\t" << crossSeveringTop << "\r\n";
    of << "crossSeveringStartAngle\t" << crossSeveringStartAngle << "\r\n";
    of << "########## Edge catastrophe specification ####\r\n";
    map<string, double>::iterator it;
    of << "########################\r\n";
    of << "edgNumber\t" << edgNumber << "\r\n";
    for (it = edgCatMap.begin(); it != edgCatMap.end(); it++)
    {
        of << it->first << "\t" << it->second << "\r\n";
    }
    of << "########################\r\n";
    of << "faceNumber\t" << faceNumber << "\r\n";
    for (it = faceCatMap.begin(); it != faceCatMap.end(); it++)
    {
        of << it->first << "\t" << it->second << "\r\n";
    }
    of << "########################\r\n";
    of << "edgeCatastropheSmooth\t" << edgeCatastropheSmooth << "\r\n";
    of << "edgeCatastropheEnabled\t" << edgeCatastropheEnabled << "\r\n";
    of << "PPBkNucFraction\t" << PPBkNucFraction << "\r\n";
    of << "proportionalCatastrophes\t" << proportionalCatastrophes << "\r\n";
    of << "########## specific conditions ###########\r\n";
    of << "restrictedPool\t" << restrictedPool << "\r\n";
    of << "interactionType\t" << interactionType << "\r\n";
    of << "bundleType\t" << bundleType << "\r\n";
    of << "########## Random seed and run time #####\r\n";
    of << "memoryLimit\t" << memoryLimit << "\r\n";
    of << "wallClockLimit\t" << wallClockLimit << "\r\n";
    of << "stopTime\t" << stopTime << "\r\n";
    of << "random_seed\t" << seed << "\r\n";
    of << "########## Movie select ############\r\n";
    of << "movieEnabled\t" << movieEnabled << "\r\n";
    of << "movieFrameInterval\t" << movieFrameInterval << "\r\n";
    of << "measurementInterval\t" << measurementInterval << "\r\n";
    of << "########## Display output ##########\r\n";
    of << "showMesh\t" << showMesh << "\r\n";
    of << "showOutput\t" << showOutput << "\r\n";
    of << "########## Geometry ################\r\n";
    of << "geometryScaleFactor\t" << geometryScaleFactor << "\r\n";
    of << "geometry\t" << geometry << "\r\n";

    return true;
}

template <class T>
inline void load_param(bool& recognized, ifstream& file, string& id, const char* tag, T& target)
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

template <class T>
inline void load_param_maps(
bool& recognized, ifstream& file, string& id, const char* tag, T& target, map<string, T>& edg_cat_map)
{
    if (!recognized)
    {
        if (id == tag)
        {
            file >> target;
            string key{ tag };
            edg_cat_map[key] = target;
            recognized = true;
        }
    }
    return;
}

template <class T>
inline void load_param_random(bool& recognized, ifstream& file, string& id, const char* tag, T& target)
{
    std::stringstream ss;
    if (!recognized)
    {
        if (id == tag)
        {
            if (id == "random_seed")
            {
                string temp;
                file >> temp;
                std::string ulong_max{ std::to_string(ULONG_MAX) };
                while (temp > ulong_max)
                {
                    cout << temp << " Bigger than" << ulong_max << endl;
                    temp = temp.substr(1);
                }
                ss.str(temp);
                ss >> target;
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

template <class T>
inline void load_param_enum_text(
bool& recognized, ifstream& file, string& id, const char* tag, T& target, string* compare_list, T compare_count)
{
    if (!recognized)
    {
        if (id == tag)
        {
            int i = 0;
            string temp;
            file >> temp;
            while (i < compare_count)
            {
                if (temp == compare_list[i])
                {
                    target = (T)i;
                    break;
                }
                i++;
            }
            if (i == compare_count)
            {
                cerr << "ERROR: No valid value for " << tag << ": " << temp << "\n";
                exit(-1);
            }
            recognized = true;
        }
    }
    return;
}

bool Parameters::read_from_file(const char* pf, bool initial_run)
{
    if ((pf == nullptr) || (*pf == 0))
    {
        return false;
    }

    ifstream par_file;
    par_file.open(pf);
    if (!par_file.good())
    {
        cerr << "Could not open parameter file [" << pf << "].\n";
        return false;
    }

    if (showOutput == 1)
    {
        cout << "Reading parameter file [" << pf << "]\n";
    }

    string id;
    bool recognized;
    while (par_file.good())
    {
        par_file >> id;
        if (!par_file.good())
        {
            break;
        }

        if (id[0] == '#')
        {
            getline(par_file, id);
        }
        else
        {
            recognized = false;

            if (initial_run)
            {
                load_param_random(recognized, par_file, id, "random_seed", seed);
            }

            load_param(recognized, par_file, id, "outputDir", outputDir);
            load_param(recognized, par_file, id, "inputDir", inputDir);
            load_param(recognized, par_file, id, "createSubdir", createSubdir);
            load_param(recognized, par_file, id, "newParameterReadInterval", newParameterReadInterval);
            load_param(recognized, par_file, id, "newParameterFile", newParameterFile);
            load_param(recognized, par_file, id, "zippering", zipperingEnabled);
            load_param(recognized, par_file, id, "ind_cat", catastrophesEnabled);
            load_param(recognized, par_file, id, "proportionalCatastrophes", proportionalCatastrophes);
            load_param(recognized, par_file, id, "induced_cat_fraction", inducedCatastropheFraction);
            load_param(recognized, par_file, id, "cat_start_angle", catStartAngle);
            load_param(recognized, par_file, id, "zipFraction", zipFraction);
            load_param(recognized, par_file, id, "magic_angle", magicAngle);
            load_param(recognized, par_file, id, "vPlus", vPlus);
            load_param(recognized, par_file, id, "vMin", vMin);
            load_param(recognized, par_file, id, "vTM", vTM);
            load_param(recognized, par_file, id, "kCat", kCat);
            load_param(recognized, par_file, id, "kRes", kRes);
            load_param(recognized, par_file, id, "kSev", kSev);
            load_param(recognized, par_file, id, "kCross", kCross);
            load_param(recognized, par_file, id, "poolDensity", poolDensity);
            load_param(recognized, par_file, id, "treadmilling", treadmillingEnabled);
            load_param(recognized, par_file, id, "severing", severingEnabled);
            load_param(recognized, par_file, id, "crossSevering", crossSeveringEnabled);
            load_param(recognized, par_file, id, "crossSeveringTop", crossSeveringTop);
            load_param(recognized, par_file, id, "crossSeveringStartAngle", crossSeveringStartAngle);
            load_param(recognized, par_file, id, "restrictedPool", restrictedPool);
            load_param(recognized, par_file, id, "edgeCatastropheEnabled", edgeCatastropheEnabled);
            load_param(recognized, par_file, id, "edgeCatastropheSmooth", edgeCatastropheSmooth);
            load_param(recognized, par_file, id, "pCatSpecialEdgeMax", pCatSpecialEdgeMax);
            load_param(recognized, par_file, id, "pCatRegularEdgeMax", pCatRegularEdgeMax);
            load_param(recognized, par_file, id, "PPB", PPB);
            load_param(recognized, par_file, id, "PPBformingTime", PPBformingTime);
            load_param(recognized, par_file, id, "kNuc", kNuc);
            load_param(recognized, par_file, id, "stopTime", stopTime);
            load_param(recognized, par_file, id, "measurementInterval", measurementInterval);
            load_param(recognized, par_file, id, "wallClockLimit", wallClockLimit);
            load_param(recognized, par_file, id, "memoryLimit", memoryLimit);
            load_param(recognized, par_file, id, "movieEnabled", movieEnabled);
            load_param(recognized, par_file, id, "movieFrameInterval", movieFrameInterval);
            load_param(recognized, par_file, id, "geometryScaleFactor", geometryScaleFactor);
            load_param_enum_text<NucleationType>(
            recognized, par_file, id, "nucleationType", nucleationType, nucleation_type_text, NUC_COUNT_LAST);
            load_param_enum_text<InteractionType>(
            recognized, par_file, id, "interactionType", interactionType, interaction_type_text, INT_COUNT_LAST);
            load_param_enum_text<BundleType>(
            recognized, par_file, id, "bundleType", bundleType, bundle_type_text, BDL_COUNT_LAST);
            load_param(recognized, par_file, id, "showOutput", showOutput);
            load_param(recognized, par_file, id, "showMesh", showMesh);
            load_param(recognized, par_file, id, "geometry", geometry);
            load_param(recognized, par_file, id, "edgNumber", edgNumber);

            for (int edg = 1; edg <= edgNumber; edg++)
            {
                double edg_cat_val;
                string edg_cat_type1 = "edgCat_" + std::to_string(edg);
                char* edg_cat_type = new char[edg_cat_type1.length() + 1];
                strcpy(edg_cat_type, edg_cat_type1.c_str());
                load_param_maps(recognized, par_file, id, edg_cat_type, edg_cat_val, edgCatMap);
                delete[] edg_cat_type;
            }

            load_param(recognized, par_file, id, "faceNumber", faceNumber);
            for (int face = 1; face <= faceNumber; face++)
            {
                double face_cat_val;
                string face_cat_type1 = "faceCat_" + std::to_string(face);
                char* face_cat_type = new char[face_cat_type1.length() + 1];
                strcpy(face_cat_type, face_cat_type1.c_str());
                load_param_maps(recognized, par_file, id, face_cat_type, face_cat_val, faceCatMap);
                delete[] face_cat_type;
            }
            if (!recognized)
            {
                cerr << "Parameter not recognized: " << id << "\n";
            }
            getline(par_file, id);
        }
    }

    if (showOutput == 1)
    {
        cout << "Finished loading parameters.\n";
    }

    return true;
}

void Parameters::verify_parameters()
{
    // adjust the parameter values in accordance with the swith (on/off) values
    vTM *= treadmillingEnabled;
    kSev *= severingEnabled;
    kCross *= crossSeveringEnabled;

    // treadmilling speed is not allowed to exceed the plus-end velocity
    if (vTM >= vPlus)
    {
        cerr
        << "ERROR: treadmilling velocity larger than growth velocity. Microtubule growth is not possible. Exiting.\n";
        exit(-1);
    }

    // minus-end velocity is not allowed to exceed the plus-end velocity
    if (vMin >= vPlus)
    {
        cerr << "ERROR: Shrinking microtubules overtake growing microtubules. This will lead to queueing problems. "
                "Exiting.\n";
        exit(-1);
    }

    // only minus-end velocity is allowed to be negative
    if ((vPlus <= 0) || (vTM < 0) || (kCat < 0) || (kRes < 0) || (kNuc < 0) || (kSev < 0) || (kCross < 0)
        || (poolDensity < 0))
    {
        cerr << "Only vMin is allowed to be negative. Exiting to avoid nonsensical results.\n";
        exit(-1);
    }

    // if there is a finite tubulin pool, check the validity of poolDensity
    if (restrictedPool != 0)
    {
        if ((poolDensity <= 0) || (poolDensity > VERY_LARGE))
        {
            cerr << "ERROR [finite tubulin pool]: maximum density is out of range (0 : VERY_LARGE) :" << poolDensity
                 << "\n";
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
        case NUC_ISOTROPIC:
            break;

        default:
            cout << "Unknown nucleationType. Aborting.\n";
            exit(-1);
            break;
    }

    switch (interactionType)
    {
        // catastrophe first, need no further conditional check
        case INT_CAT_FIRST:
            break;

        // zip first needs further check
        case INT_ZIP_FIRST:

            // zippering fraction must be in the range: (0.0 - 1.0)
            if (zipperingEnabled)
            {
                if (zipFraction < 0.)
                {
                    zipFraction = 0.;
                }
                if (zipFraction > 1.)
                {
                    zipFraction = 1.;
                }
            }
            break;

        case INT_COUNT_LAST:
            break;
    }

    return;
}

bool Parameters::calc_theory_parameters(void)
{
    // convert the angle values from degrees to radians
    double ma_rad = magicAngle * PI / 180.0;
    double cs_rad = catStartAngle * PI / 180.0;

    switch (interactionType)
    {
        // zippering first is requested
        case INT_ZIP_FIRST:

            // fzippering is enabled
            if (zipperingEnabled)
            {
                z0calc = (4. / PI) * (1. - cos(ma_rad)) * zipFraction;
            }

            // zippering is disabled
            else
            {
                z0calc = 0;
            }

            // catastrophe is enabled
            if (catastrophesEnabled)
            {
                // proportional catastrophe is enabled
                if (proportionalCatastrophes)
                {
                    // magic angle > catastrophe start angle
                    if (magicAngle > catStartAngle)
                    {
                        c0calc = (4.0 / PI) * 2.0 * inducedCatastropheFraction
                                 * (1. - sin(ma_rad) + (ma_rad - cs_rad) * cos(ma_rad)) / (PI - 2. * cs_rad);
                    }

                    // magic angle <= catastrophe start angle
                    else
                    {
                        c0calc
                        = (4.0 / PI) * 2.0 * inducedCatastropheFraction * (1.0 - sin(cs_rad)) / (PI - 2.0 * cs_rad);
                    }
                }

                // proportional catastrophe is disabled
                else
                {
                    c0calc = (4.0 / PI) * inducedCatastropheFraction * cos(max(ma_rad, cs_rad));
                }
            }

            // catastrophe is disabled
            else
            {
                c0calc = 0;
            }

            // calculate x0, using the values of c0 and/or z0
            x0calc = (4.0 / PI) - z0calc - c0calc;
            break;

        // catastrophe first is requested
        case INT_CAT_FIRST:

            // catastrophe is enabled
            if (catastrophesEnabled)
            {
                // proportional catastrophe is enabled
                if (proportionalCatastrophes)
                {
                    double temp = max(0., cs_rad);
                    c0calc = (4.0 / PI) * 2.0 * inducedCatastropheFraction
                             * (1.0 - sin(temp) + (temp - cs_rad) * cos(temp)) / (PI - 2.0 * cs_rad);
                }

                // proportional catastrophe is disabled
                else
                {
                    c0calc = (4.0 / PI) * inducedCatastropheFraction * cos(max(cs_rad, 0.));
                }
            }

            // zippering is enabled
            if (zipperingEnabled)
            {
                z0calc = (4.0 / PI) * (1.0 - cos(ma_rad));

                // catastrophe is enabled and  (magic angle > catastrophe start angle)
                if (catastrophesEnabled && (cs_rad < ma_rad))
                {
                    // proportional catastrophe is enabled
                    if (proportionalCatastrophes)
                    {
                        double left_boundary = max(0., cs_rad);
                        z0calc -= (4.0 / PI) * 2.0 * inducedCatastropheFraction
                                  * ((left_boundary - cs_rad) * cos(left_boundary) - (ma_rad - cs_rad) * cos(ma_rad)
                                     - sin(left_boundary) + sin(ma_rad))
                                  / (PI - 2.0 * cs_rad);
                    }

                    // proportional catastrophe is disabled
                    else
                    {
                        z0calc -= (4.0 / PI) * inducedCatastropheFraction * (cos(cs_rad) - cos(ma_rad));
                    }
                }
            }

            // zippering is disabled
            else
            {
                z0calc = 0;
            }

            // calculate x0, using the values of c0 and/or z0
            x0calc = (4.0 / PI) - z0calc - c0calc;
            break;

        default:
            cerr << "no function to calculate z0, c0, x0 etc. (interactionType = " << interactionType << ")\n";
            return false;
    }

    // theoretial value of G
    double l0 = pow((2.0 * (-vMin + vTM) * (vPlus - vTM) * (vPlus - vTM)) / (kNuc * vPlus * (-vMin + vPlus)), 1. / 3.);
    double lavg_inv = kRes / (-vMin + vTM) - kCat / (vPlus - vTM);
    double g_0 = l0 * lavg_inv;

    if (showOutput == 1)
    {
        cout << "===============================================" << endl;
        cout << " Calculated theory parameters" << endl;
        cout << "===============================================" << endl;
        cout << "l0 = " << l0 << endl;
        ;
        cout << "lavg = " << -1. / lavg_inv << endl;
        cout << "G_0 = -l0/lavg = " << g_0 << endl;
        cout << "===============================================" << endl;
    }

    return (true);
}
