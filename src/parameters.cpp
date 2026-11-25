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
    v_plus(0.08),
    v_min(-0.16),
    v_tm(0.01),
    k_sev(0.001),
    k_cross(0.001),
    k_cat(0.005),
    k_res(0.007),
    k_nuc(0.001),
    pp_bforming_time(VERY_LARGE),
    pool_density(5.0),
    treadmilling_enabled(1),
    severing_enabled(0),
    cross_severing_enabled(0),
    cross_severing_top(1),
    cross_severing_start_angle(0),
    ppb(0),
    restricted_pool(0),
    edge_catastrophe_enabled(0),
    edge_catastrophe_smooth(0),
    p_cat_special_edge_max(0.0),
    p_cat_regular_edge_max(0.0),
    edg_number(0),
    face_number(0),
    nucleation_type(NucleationType::ISOTROPIC),
    zippering_enabled(1),
    catastrophes_enabled(1),
    proportional_catastrophes(0),
    induced_catastrophe_fraction(0.5),
    zip_fraction(1.0),
    pp_bk_nuc_fraction(0.0),
    cat_start_angle(0),
    magic_angle(40),
    interaction_type(InteractionType::ZIP_FIRST),
    bundle_type(BundleType::SIMPLE),
    geom_param("0"),
    seed(time(nullptr)),
    stop_time(36000),
    measurement_interval(600),
    wall_clock_limit(1000),
    memory_limit(10000000),
    input_dir("."),
    output_dir("."),
    create_subdir(1),
    new_parameter_read_interval(VERY_LARGE),
    new_parameter_file(""),
    movie_enabled(0),
    show_output(0),
    show_mesh(0),
    geometry("2D-plane_1_0_0"),
    movie_frame_interval(-1),
    geometry_scale_factor(15.0),
    c0_calc(0),
    x0_calc(0),
    z0_calc(0)
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
    ofstream& of = system->parameter_file;
    of << "# CorticalSim parameters\r\n";
    of << "########## Input-Output path ########################\r\n";
    of << "inputDir\t" << input_dir << "\r\n";
    of << "outputDir\t" << output_dir << "\r\n";
    of << "newParameterReadInterval\t" << new_parameter_read_interval << "\r\n";
    of << "createSubdir\t" << create_subdir << "\r\n";
    of << "newParameterFile\t" << new_parameter_file << "\r\n";
    of << "########## Dynamic parameters ##########\r\n";
    of << "kSev\t" << k_sev << "\r\n";
    of << "kCross\t" << k_cross << "\r\n";
    of << "vPlus\t" << v_plus << "\r\n";
    of << "vMin\t" << v_min << "\r\n";
    of << "vTM\t" << v_tm << "\r\n";
    of << "kCat\t" << k_cat << "\r\n";
    of << "kRes\t" << k_res << "\r\n";
    of << "kNuc\t" << k_nuc << "\r\n";
    of << "magic_angle\t" << magic_angle << "\r\n";
    of << "poolDensity\t" << pool_density << "\r\n";
    of << "nucleationType\t" << static_cast<int>(nucleation_type) << "\r\n";
    of << "########## Zip-IndCat-Severing events ######\r\n";
    of << "treadmilling\t" << treadmilling_enabled << "\r\n";
    of << "ind_cat\t" << catastrophes_enabled << "\r\n";
    of << "cat_start_angle\t" << cat_start_angle << "\r\n";
    of << "induced_cat_fraction\t" << induced_catastrophe_fraction << "\r\n";
    of << "zippering\t" << zippering_enabled << "\r\n";
    of << "zipFraction\t" << zip_fraction << "\r\n";
    of << "severing\t" << severing_enabled << "\r\n";
    of << "crossSevering\t" << cross_severing_enabled << "\r\n";
    of << "crossSeveringTop\t" << cross_severing_top << "\r\n";
    of << "crossSeveringStartAngle\t" << cross_severing_start_angle << "\r\n";
    of << "########## Edge catastrophe specification ####\r\n";
    map<string, double>::iterator it;
    of << "########################\r\n";
    of << "edgNumber\t" << edg_number << "\r\n";
    for (it = edg_cat_map.begin(); it != edg_cat_map.end(); it++)
    {
        of << it->first << "\t" << it->second << "\r\n";
    }
    of << "########################\r\n";
    of << "faceNumber\t" << face_number << "\r\n";
    for (it = face_cat_map.begin(); it != face_cat_map.end(); it++)
    {
        of << it->first << "\t" << it->second << "\r\n";
    }
    of << "########################\r\n";
    of << "edgeCatastropheSmooth\t" << edge_catastrophe_smooth << "\r\n";
    of << "edgeCatastropheEnabled\t" << edge_catastrophe_enabled << "\r\n";
    of << "PPBkNucFraction\t" << pp_bk_nuc_fraction << "\r\n";
    of << "proportionalCatastrophes\t" << proportional_catastrophes << "\r\n";
    of << "########## specific conditions ###########\r\n";
    of << "restrictedPool\t" << restricted_pool << "\r\n";
    of << "interactionType\t" << static_cast<int>(interaction_type) << "\r\n";
    of << "bundleType\t" << static_cast<int>(bundle_type) << "\r\n";
    of << "########## Random seed and run time #####\r\n";
    of << "memoryLimit\t" << memory_limit << "\r\n";
    of << "wallClockLimit\t" << wall_clock_limit << "\r\n";
    of << "stopTime\t" << stop_time << "\r\n";
    of << "random_seed\t" << seed << "\r\n";
    of << "########## Movie select ############\r\n";
    of << "movieEnabled\t" << movie_enabled << "\r\n";
    of << "movieFrameInterval\t" << movie_frame_interval << "\r\n";
    of << "measurementInterval\t" << measurement_interval << "\r\n";
    of << "########## Display output ##########\r\n";
    of << "showMesh\t" << show_mesh << "\r\n";
    of << "showOutput\t" << show_output << "\r\n";
    of << "########## Geometry ################\r\n";
    of << "geometryScaleFactor\t" << geometry_scale_factor << "\r\n";
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
            // TODO: Fragile - refactor as a priority.
            while (i < static_cast<int>(compare_count))
            {
                if (temp == compare_list[i])
                {
                    target = (T)i;
                    break;
                }
                i++;
            }
            // TODO: Fragile - refactor as a priority.
            if (i == static_cast<int>(compare_count))
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

    if (show_output == 1)
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

            load_param(recognized, par_file, id, "outputDir", output_dir);
            load_param(recognized, par_file, id, "inputDir", input_dir);
            load_param(recognized, par_file, id, "createSubdir", create_subdir);
            load_param(recognized, par_file, id, "newParameterReadInterval", new_parameter_read_interval);
            load_param(recognized, par_file, id, "newParameterFile", new_parameter_file);
            load_param(recognized, par_file, id, "zippering", zippering_enabled);
            load_param(recognized, par_file, id, "ind_cat", catastrophes_enabled);
            load_param(recognized, par_file, id, "proportionalCatastrophes", proportional_catastrophes);
            load_param(recognized, par_file, id, "induced_cat_fraction", induced_catastrophe_fraction);
            load_param(recognized, par_file, id, "cat_start_angle", cat_start_angle);
            load_param(recognized, par_file, id, "zipFraction", zip_fraction);
            load_param(recognized, par_file, id, "magic_angle", magic_angle);
            load_param(recognized, par_file, id, "vPlus", v_plus);
            load_param(recognized, par_file, id, "vMin", v_min);
            load_param(recognized, par_file, id, "vTM", v_tm);
            load_param(recognized, par_file, id, "kCat", k_cat);
            load_param(recognized, par_file, id, "kRes", k_res);
            load_param(recognized, par_file, id, "kSev", k_sev);
            load_param(recognized, par_file, id, "kCross", k_cross);
            load_param(recognized, par_file, id, "poolDensity", pool_density);
            load_param(recognized, par_file, id, "treadmilling", treadmilling_enabled);
            load_param(recognized, par_file, id, "severing", severing_enabled);
            load_param(recognized, par_file, id, "crossSevering", cross_severing_enabled);
            load_param(recognized, par_file, id, "crossSeveringTop", cross_severing_top);
            load_param(recognized, par_file, id, "crossSeveringStartAngle", cross_severing_start_angle);
            load_param(recognized, par_file, id, "restrictedPool", restricted_pool);
            load_param(recognized, par_file, id, "edgeCatastropheEnabled", edge_catastrophe_enabled);
            load_param(recognized, par_file, id, "edgeCatastropheSmooth", edge_catastrophe_smooth);
            load_param(recognized, par_file, id, "pCatSpecialEdgeMax", p_cat_special_edge_max);
            load_param(recognized, par_file, id, "pCatRegularEdgeMax", p_cat_regular_edge_max);
            load_param(recognized, par_file, id, "PPB", ppb);
            load_param(recognized, par_file, id, "PPBformingTime", pp_bforming_time);
            load_param(recognized, par_file, id, "kNuc", k_nuc);
            load_param(recognized, par_file, id, "stopTime", stop_time);
            load_param(recognized, par_file, id, "measurementInterval", measurement_interval);
            load_param(recognized, par_file, id, "wallClockLimit", wall_clock_limit);
            load_param(recognized, par_file, id, "memoryLimit", memory_limit);
            load_param(recognized, par_file, id, "movieEnabled", movie_enabled);
            load_param(recognized, par_file, id, "movieFrameInterval", movie_frame_interval);
            load_param(recognized, par_file, id, "geometryScaleFactor", geometry_scale_factor);
            load_param_enum_text<NucleationType>(recognized,
                                                 par_file,
                                                 id,
                                                 "nucleationType",
                                                 nucleation_type,
                                                 nucleation_type_text,
                                                 NucleationType::COUNT_LAST);
            load_param_enum_text<InteractionType>(recognized,
                                                  par_file,
                                                  id,
                                                  "interactionType",
                                                  interaction_type,
                                                  interaction_type_text,
                                                  InteractionType::COUNT_LAST);
            load_param_enum_text<BundleType>(
            recognized, par_file, id, "bundleType", bundle_type, bundle_type_text, BundleType::COUNT_LAST);
            load_param(recognized, par_file, id, "showOutput", show_output);
            load_param(recognized, par_file, id, "showMesh", show_mesh);
            load_param(recognized, par_file, id, "geometry", geometry);
            load_param(recognized, par_file, id, "edgNumber", edg_number);

            for (int edg = 1; edg <= edg_number; edg++)
            {
                double edg_cat_val;
                string edg_cat_type1 = "edgCat_" + std::to_string(edg);
                char* edg_cat_type = new char[edg_cat_type1.length() + 1];
                strcpy(edg_cat_type, edg_cat_type1.c_str());
                load_param_maps(recognized, par_file, id, edg_cat_type, edg_cat_val, edg_cat_map);
                delete[] edg_cat_type;
            }

            load_param(recognized, par_file, id, "faceNumber", face_number);
            for (int face = 1; face <= face_number; face++)
            {
                double face_cat_val;
                string face_cat_type1 = "faceCat_" + std::to_string(face);
                char* face_cat_type = new char[face_cat_type1.length() + 1];
                strcpy(face_cat_type, face_cat_type1.c_str());
                load_param_maps(recognized, par_file, id, face_cat_type, face_cat_val, face_cat_map);
                delete[] face_cat_type;
            }
            if (!recognized)
            {
                cerr << "Parameter not recognized: " << id << "\n";
            }
            getline(par_file, id);
        }
    }

    if (show_output == 1)
    {
        cout << "Finished loading parameters.\n";
    }

    return true;
}

void Parameters::verify_parameters()
{
    // adjust the parameter values in accordance with the swith (on/off) values
    v_tm *= treadmilling_enabled;
    k_sev *= severing_enabled;
    k_cross *= cross_severing_enabled;

    // treadmilling speed is not allowed to exceed the plus-end velocity
    if (v_tm >= v_plus)
    {
        cerr
        << "ERROR: treadmilling velocity larger than growth velocity. Microtubule growth is not possible. Exiting.\n";
        exit(-1);
    }

    // minus-end velocity is not allowed to exceed the plus-end velocity
    if (v_min >= v_plus)
    {
        cerr << "ERROR: Shrinking microtubules overtake growing microtubules. This will lead to queueing problems. "
                "Exiting.\n";
        exit(-1);
    }

    // only minus-end velocity is allowed to be negative
    if ((v_plus <= 0) || (v_tm < 0) || (k_cat < 0) || (k_res < 0) || (k_nuc < 0) || (k_sev < 0) || (k_cross < 0)
        || (pool_density < 0))
    {
        cerr << "Only vMin is allowed to be negative. Exiting to avoid nonsensical results.\n";
        exit(-1);
    }

    // if there is a finite tubulin pool, check the validity of poolDensity
    if (restricted_pool != 0)
    {
        if ((pool_density <= 0) || (pool_density > VERY_LARGE))
        {
            cerr << "ERROR [finite tubulin pool]: maximum density is out of range (0 : VERY_LARGE) :" << pool_density
                 << "\n";
            exit(-1);
        }
    }

    // it is not possible to create movie without valid movie-interval
    if (movie_enabled)
    {
        if (movie_frame_interval < ZERO_CUTOFF)
        {
            cerr << "ERROR: movie recording cannot be enabled without specifying 'movieFrameInterval'\n";
            exit(-1);
        }
    }

    switch (nucleation_type)
    {
            // isotropic nucleation need no conditional check
        case NucleationType::ISOTROPIC:
            break;

        default:
            cout << "Unknown nucleationType. Aborting.\n";
            exit(-1);
            break;
    }

    switch (interaction_type)
    {
        // catastrophe first, need no further conditional check
        case InteractionType::CAT_FIRST:
            break;

        // zip first needs further check
        case InteractionType::ZIP_FIRST:

            // zippering fraction must be in the range: (0.0 - 1.0)
            if (zippering_enabled)
            {
                if (zip_fraction < 0.)
                {
                    zip_fraction = 0.;
                }
                if (zip_fraction > 1.)
                {
                    zip_fraction = 1.;
                }
            }
            break;

        case InteractionType::COUNT_LAST:
            break;
    }

    return;
}

bool Parameters::calc_theory_parameters(void)
{
    // convert the angle values from degrees to radians
    double ma_rad = magic_angle * PI / 180.0;
    double cs_rad = cat_start_angle * PI / 180.0;

    switch (interaction_type)
    {
        // zippering first is requested
        case InteractionType::ZIP_FIRST:

            // fzippering is enabled
            if (zippering_enabled)
            {
                z0_calc = (4. / PI) * (1. - cos(ma_rad)) * zip_fraction;
            }

            // zippering is disabled
            else
            {
                z0_calc = 0;
            }

            // catastrophe is enabled
            if (catastrophes_enabled)
            {
                // proportional catastrophe is enabled
                if (proportional_catastrophes)
                {
                    // magic angle > catastrophe start angle
                    if (magic_angle > cat_start_angle)
                    {
                        c0_calc = (4.0 / PI) * 2.0 * induced_catastrophe_fraction
                                  * (1. - sin(ma_rad) + (ma_rad - cs_rad) * cos(ma_rad)) / (PI - 2. * cs_rad);
                    }

                    // magic angle <= catastrophe start angle
                    else
                    {
                        c0_calc
                        = (4.0 / PI) * 2.0 * induced_catastrophe_fraction * (1.0 - sin(cs_rad)) / (PI - 2.0 * cs_rad);
                    }
                }

                // proportional catastrophe is disabled
                else
                {
                    c0_calc = (4.0 / PI) * induced_catastrophe_fraction * cos(max(ma_rad, cs_rad));
                }
            }

            // catastrophe is disabled
            else
            {
                c0_calc = 0;
            }

            // calculate x0, using the values of c0 and/or z0
            x0_calc = (4.0 / PI) - z0_calc - c0_calc;
            break;

        // catastrophe first is requested
        case InteractionType::CAT_FIRST:

            // catastrophe is enabled
            if (catastrophes_enabled)
            {
                // proportional catastrophe is enabled
                if (proportional_catastrophes)
                {
                    double temp = max(0., cs_rad);
                    c0_calc = (4.0 / PI) * 2.0 * induced_catastrophe_fraction
                              * (1.0 - sin(temp) + (temp - cs_rad) * cos(temp)) / (PI - 2.0 * cs_rad);
                }

                // proportional catastrophe is disabled
                else
                {
                    c0_calc = (4.0 / PI) * induced_catastrophe_fraction * cos(max(cs_rad, 0.));
                }
            }

            // zippering is enabled
            if (zippering_enabled)
            {
                z0_calc = (4.0 / PI) * (1.0 - cos(ma_rad));

                // catastrophe is enabled and  (magic angle > catastrophe start angle)
                if (catastrophes_enabled && (cs_rad < ma_rad))
                {
                    // proportional catastrophe is enabled
                    if (proportional_catastrophes)
                    {
                        double left_boundary = max(0., cs_rad);
                        z0_calc -= (4.0 / PI) * 2.0 * induced_catastrophe_fraction
                                   * ((left_boundary - cs_rad) * cos(left_boundary) - (ma_rad - cs_rad) * cos(ma_rad)
                                      - sin(left_boundary) + sin(ma_rad))
                                   / (PI - 2.0 * cs_rad);
                    }

                    // proportional catastrophe is disabled
                    else
                    {
                        z0_calc -= (4.0 / PI) * induced_catastrophe_fraction * (cos(cs_rad) - cos(ma_rad));
                    }
                }
            }

            // zippering is disabled
            else
            {
                z0_calc = 0;
            }

            // calculate x0, using the values of c0 and/or z0
            x0_calc = (4.0 / PI) - z0_calc - c0_calc;
            break;

        default:
            cerr << "no function to calculate z0, c0, x0 etc. (interactionType = " << static_cast<int>(interaction_type)
                 << ")\n";
            return false;
    }

    // theoretial value of G
    double l0
    = pow((2.0 * (-v_min + v_tm) * (v_plus - v_tm) * (v_plus - v_tm)) / (k_nuc * v_plus * (-v_min + v_plus)), 1. / 3.);
    double lavg_inv = k_res / (-v_min + v_tm) - k_cat / (v_plus - v_tm);
    double g_0 = l0 * lavg_inv;

    if (show_output == 1)
    {
        cout << "===============================================" << endl;
        cout << " Calculated theory parameters" << endl;
        cout << "===============================================" << endl;
        cout << "l0 = " << l0 << endl;
        cout << "lavg = " << -1. / lavg_inv << endl;
        cout << "G_0 = -l0/lavg = " << g_0 << endl;
        cout << "===============================================" << endl;
    }

    return (true);
}
