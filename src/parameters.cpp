#include <sstream>

#include "parameters.hpp"
#include "system.hpp"

string direction_type_text[] = { "forward", "backward" };
string geometry_type_text[] = { "cell" };
string interaction_type_text[] = { "zipFirst", "catFirst", "minimalFourier" };
string nucleation_type_text[] = { "isotropic", "biased", "ellipse" };
string bundle_type_text[] = { "simple", "sticky", "noZip", "multiCollision", "Ncollision" };

Parameters::Parameters(System* s):
    m_system(s),
    m_v_plus(0.08),
    m_v_min(-0.16),
    m_v_tm(0.01),
    m_k_sev(0.001),
    m_k_cross(0.001),
    m_k_cat(0.005),
    m_k_res(0.007),
    m_k_nuc(0.001),
    m_pp_bforming_time(VERY_LARGE),
    m_pool_density(5.0),
    m_treadmilling_enabled(1),
    m_severing_enabled(0),
    m_cross_severing_enabled(0),
    m_cross_severing_top(1),
    m_cross_severing_start_angle(0),
    m_ppb(0),
    m_restricted_pool(0),
    m_edge_catastrophe_enabled(0),
    m_edge_catastrophe_smooth(0),
    m_p_cat_special_edge_max(0.0),
    m_p_cat_regular_edge_max(0.0),
    m_edg_number(0),
    m_face_number(0),
    m_nucleation_type(NUC_ISOTROPIC),
    m_zippering_enabled(1),
    m_catastrophes_enabled(1),
    m_proportional_catastrophes(0),
    m_induced_catastrophe_fraction(0.5),
    m_zip_fraction(1.0),
    m_pp_bk_nuc_fraction(0.0),
    m_cat_start_angle(0),
    m_magic_angle(40),
    m_interaction_type(INT_ZIP_FIRST),
    m_bundle_type(BDL_SIMPLE),
    m_geom_param("0"),
    m_seed(time(nullptr)),
    m_stop_time(36000),
    m_measurement_interval(600),
    m_wall_clock_limit(1000),
    m_memory_limit(10000000),
    m_input_dir("."),
    m_output_dir("."),
    m_create_subdir(1),
    m_new_parameter_read_interval(VERY_LARGE),
    m_new_parameter_file(""),
    m_movie_enabled(0),
    m_show_output(0),
    m_show_mesh(0),
    m_geometry("2D-plane_1_0_0"),
    m_movie_frame_interval(-1),
    m_geometry_scale_factor(15.0),
    m_c0calc(0),
    m_x0calc(0),
    m_z0calc(0)
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
    ofstream& of = m_system->m_parameter_file;
    of << "# CorticalSim parameters\r\n";
    of << "########## Input-Output path ########################\r\n";
    of << "inputDir\t" << m_input_dir << "\r\n";
    of << "outputDir\t" << m_output_dir << "\r\n";
    of << "newParameterReadInterval\t" << m_new_parameter_read_interval << "\r\n";
    of << "createSubdir\t" << m_create_subdir << "\r\n";
    of << "newParameterFile\t" << m_new_parameter_file << "\r\n";
    of << "########## Dynamic parameters ##########\r\n";
    of << "kSev\t" << m_k_sev << "\r\n";
    of << "kCross\t" << m_k_cross << "\r\n";
    of << "vPlus\t" << m_v_plus << "\r\n";
    of << "vMin\t" << m_v_min << "\r\n";
    of << "vTM\t" << m_v_tm << "\r\n";
    of << "kCat\t" << m_k_cat << "\r\n";
    of << "kRes\t" << m_k_res << "\r\n";
    of << "kNuc\t" << m_k_nuc << "\r\n";
    of << "magic_angle\t" << m_magic_angle << "\r\n";
    of << "poolDensity\t" << m_pool_density << "\r\n";
    of << "nucleationType\t" << m_nucleation_type << "\r\n";
    of << "########## Zip-IndCat-Severing events ######\r\n";
    of << "treadmilling\t" << m_treadmilling_enabled << "\r\n";
    of << "ind_cat\t" << m_catastrophes_enabled << "\r\n";
    of << "cat_start_angle\t" << m_cat_start_angle << "\r\n";
    of << "induced_cat_fraction\t" << m_induced_catastrophe_fraction << "\r\n";
    of << "zippering\t" << m_zippering_enabled << "\r\n";
    of << "zipFraction\t" << m_zip_fraction << "\r\n";
    of << "severing\t" << m_severing_enabled << "\r\n";
    of << "crossSevering\t" << m_cross_severing_enabled << "\r\n";
    of << "crossSeveringTop\t" << m_cross_severing_top << "\r\n";
    of << "crossSeveringStartAngle\t" << m_cross_severing_start_angle << "\r\n";
    of << "########## Edge catastrophe specification ####\r\n";
    map<string, double>::iterator it;
    of << "########################\r\n";
    of << "edgNumber\t" << m_edg_number << "\r\n";
    for (it = m_edg_cat_map.begin(); it != m_edg_cat_map.end(); it++)
    {
        of << it->first << "\t" << it->second << "\r\n";
    }
    of << "########################\r\n";
    of << "faceNumber\t" << m_face_number << "\r\n";
    for (it = m_face_cat_map.begin(); it != m_face_cat_map.end(); it++)
    {
        of << it->first << "\t" << it->second << "\r\n";
    }
    of << "########################\r\n";
    of << "edgeCatastropheSmooth\t" << m_edge_catastrophe_smooth << "\r\n";
    of << "edgeCatastropheEnabled\t" << m_edge_catastrophe_enabled << "\r\n";
    of << "PPBkNucFraction\t" << m_pp_bk_nuc_fraction << "\r\n";
    of << "proportionalCatastrophes\t" << m_proportional_catastrophes << "\r\n";
    of << "########## specific conditions ###########\r\n";
    of << "restrictedPool\t" << m_restricted_pool << "\r\n";
    of << "interactionType\t" << m_interaction_type << "\r\n";
    of << "bundleType\t" << m_bundle_type << "\r\n";
    of << "########## Random seed and run time #####\r\n";
    of << "memoryLimit\t" << m_memory_limit << "\r\n";
    of << "wallClockLimit\t" << m_wall_clock_limit << "\r\n";
    of << "stopTime\t" << m_stop_time << "\r\n";
    of << "random_seed\t" << m_seed << "\r\n";
    of << "########## Movie select ############\r\n";
    of << "movieEnabled\t" << m_movie_enabled << "\r\n";
    of << "movieFrameInterval\t" << m_movie_frame_interval << "\r\n";
    of << "measurementInterval\t" << m_measurement_interval << "\r\n";
    of << "########## Display output ##########\r\n";
    of << "showMesh\t" << m_show_mesh << "\r\n";
    of << "showOutput\t" << m_show_output << "\r\n";
    of << "########## Geometry ################\r\n";
    of << "geometryScaleFactor\t" << m_geometry_scale_factor << "\r\n";
    of << "geometry\t" << m_geometry << "\r\n";

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

    if (m_show_output == 1)
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
                load_param_random(recognized, par_file, id, "random_seed", m_seed);
            }

            load_param(recognized, par_file, id, "outputDir", m_output_dir);
            load_param(recognized, par_file, id, "inputDir", m_input_dir);
            load_param(recognized, par_file, id, "createSubdir", m_create_subdir);
            load_param(recognized, par_file, id, "newParameterReadInterval", m_new_parameter_read_interval);
            load_param(recognized, par_file, id, "newParameterFile", m_new_parameter_file);
            load_param(recognized, par_file, id, "zippering", m_zippering_enabled);
            load_param(recognized, par_file, id, "ind_cat", m_catastrophes_enabled);
            load_param(recognized, par_file, id, "proportionalCatastrophes", m_proportional_catastrophes);
            load_param(recognized, par_file, id, "induced_cat_fraction", m_induced_catastrophe_fraction);
            load_param(recognized, par_file, id, "cat_start_angle", m_cat_start_angle);
            load_param(recognized, par_file, id, "zipFraction", m_zip_fraction);
            load_param(recognized, par_file, id, "magic_angle", m_magic_angle);
            load_param(recognized, par_file, id, "vPlus", m_v_plus);
            load_param(recognized, par_file, id, "vMin", m_v_min);
            load_param(recognized, par_file, id, "vTM", m_v_tm);
            load_param(recognized, par_file, id, "kCat", m_k_cat);
            load_param(recognized, par_file, id, "kRes", m_k_res);
            load_param(recognized, par_file, id, "kSev", m_k_sev);
            load_param(recognized, par_file, id, "kCross", m_k_cross);
            load_param(recognized, par_file, id, "poolDensity", m_pool_density);
            load_param(recognized, par_file, id, "treadmilling", m_treadmilling_enabled);
            load_param(recognized, par_file, id, "severing", m_severing_enabled);
            load_param(recognized, par_file, id, "crossSevering", m_cross_severing_enabled);
            load_param(recognized, par_file, id, "crossSeveringTop", m_cross_severing_top);
            load_param(recognized, par_file, id, "crossSeveringStartAngle", m_cross_severing_start_angle);
            load_param(recognized, par_file, id, "restrictedPool", m_restricted_pool);
            load_param(recognized, par_file, id, "edgeCatastropheEnabled", m_edge_catastrophe_enabled);
            load_param(recognized, par_file, id, "edgeCatastropheSmooth", m_edge_catastrophe_smooth);
            load_param(recognized, par_file, id, "pCatSpecialEdgeMax", m_p_cat_special_edge_max);
            load_param(recognized, par_file, id, "pCatRegularEdgeMax", m_p_cat_regular_edge_max);
            load_param(recognized, par_file, id, "PPB", m_ppb);
            load_param(recognized, par_file, id, "PPBformingTime", m_pp_bforming_time);
            load_param(recognized, par_file, id, "kNuc", m_k_nuc);
            load_param(recognized, par_file, id, "stopTime", m_stop_time);
            load_param(recognized, par_file, id, "measurementInterval", m_measurement_interval);
            load_param(recognized, par_file, id, "wallClockLimit", m_wall_clock_limit);
            load_param(recognized, par_file, id, "memoryLimit", m_memory_limit);
            load_param(recognized, par_file, id, "movieEnabled", m_movie_enabled);
            load_param(recognized, par_file, id, "movieFrameInterval", m_movie_frame_interval);
            load_param(recognized, par_file, id, "geometryScaleFactor", m_geometry_scale_factor);
            load_param_enum_text<NucleationType>(
            recognized, par_file, id, "nucleationType", m_nucleation_type, nucleation_type_text, NUC_COUNT_LAST);
            load_param_enum_text<InteractionType>(
            recognized, par_file, id, "interactionType", m_interaction_type, interaction_type_text, INT_COUNT_LAST);
            load_param_enum_text<BundleType>(
            recognized, par_file, id, "bundleType", m_bundle_type, bundle_type_text, BDL_COUNT_LAST);
            load_param(recognized, par_file, id, "showOutput", m_show_output);
            load_param(recognized, par_file, id, "showMesh", m_show_mesh);
            load_param(recognized, par_file, id, "geometry", m_geometry);
            load_param(recognized, par_file, id, "edgNumber", m_edg_number);

            for (int edg = 1; edg <= m_edg_number; edg++)
            {
                double edg_cat_val;
                string edg_cat_type1 = "edgCat_" + std::to_string(edg);
                char* edg_cat_type = new char[edg_cat_type1.length() + 1];
                strcpy(edg_cat_type, edg_cat_type1.c_str());
                load_param_maps(recognized, par_file, id, edg_cat_type, edg_cat_val, m_edg_cat_map);
                delete[] edg_cat_type;
            }

            load_param(recognized, par_file, id, "faceNumber", m_face_number);
            for (int face = 1; face <= m_face_number; face++)
            {
                double face_cat_val;
                string face_cat_type1 = "faceCat_" + std::to_string(face);
                char* face_cat_type = new char[face_cat_type1.length() + 1];
                strcpy(face_cat_type, face_cat_type1.c_str());
                load_param_maps(recognized, par_file, id, face_cat_type, face_cat_val, m_face_cat_map);
                delete[] face_cat_type;
            }
            if (!recognized)
            {
                cerr << "Parameter not recognized: " << id << "\n";
            }
            getline(par_file, id);
        }
    }

    if (m_show_output == 1)
    {
        cout << "Finished loading parameters.\n";
    }

    return true;
}

void Parameters::verify_parameters()
{
    // adjust the parameter values in accordance with the swith (on/off) values
    m_v_tm *= m_treadmilling_enabled;
    m_k_sev *= m_severing_enabled;
    m_k_cross *= m_cross_severing_enabled;

    // treadmilling speed is not allowed to exceed the plus-end velocity
    if (m_v_tm >= m_v_plus)
    {
        cerr
        << "ERROR: treadmilling velocity larger than growth velocity. Microtubule growth is not possible. Exiting.\n";
        exit(-1);
    }

    // minus-end velocity is not allowed to exceed the plus-end velocity
    if (m_v_min >= m_v_plus)
    {
        cerr << "ERROR: Shrinking microtubules overtake growing microtubules. This will lead to queueing problems. "
                "Exiting.\n";
        exit(-1);
    }

    // only minus-end velocity is allowed to be negative
    if ((m_v_plus <= 0) || (m_v_tm < 0) || (m_k_cat < 0) || (m_k_res < 0) || (m_k_nuc < 0) || (m_k_sev < 0) || (m_k_cross < 0)
        || (m_pool_density < 0))
    {
        cerr << "Only vMin is allowed to be negative. Exiting to avoid nonsensical results.\n";
        exit(-1);
    }

    // if there is a finite tubulin pool, check the validity of poolDensity
    if (m_restricted_pool != 0)
    {
        if ((m_pool_density <= 0) || (m_pool_density > VERY_LARGE))
        {
            cerr << "ERROR [finite tubulin pool]: maximum density is out of range (0 : VERY_LARGE) :" << m_pool_density
                 << "\n";
            exit(-1);
        }
    }

    // it is not possible to create movie without valid movie-interval
    if (m_movie_enabled)
    {
        if (m_movie_frame_interval < ZERO_CUTOFF)
        {
            cerr << "ERROR: movie recording cannot be enabled without specifying 'movieFrameInterval'\n";
            exit(-1);
        }
    }

    switch (m_nucleation_type)
    {
            // isotropic nucleation need no conditional check
        case NUC_ISOTROPIC:
            break;

        default:
            cout << "Unknown nucleationType. Aborting.\n";
            exit(-1);
            break;
    }

    switch (m_interaction_type)
    {
        // catastrophe first, need no further conditional check
        case INT_CAT_FIRST:
            break;

        // zip first needs further check
        case INT_ZIP_FIRST:

            // zippering fraction must be in the range: (0.0 - 1.0)
            if (m_zippering_enabled)
            {
                if (m_zip_fraction < 0.)
                {
                    m_zip_fraction = 0.;
                }
                if (m_zip_fraction > 1.)
                {
                    m_zip_fraction = 1.;
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
    double ma_rad = m_magic_angle * PI / 180.0;
    double cs_rad = m_cat_start_angle * PI / 180.0;

    switch (m_interaction_type)
    {
        // zippering first is requested
        case INT_ZIP_FIRST:

            // fzippering is enabled
            if (m_zippering_enabled)
            {
                m_z0calc = (4. / PI) * (1. - cos(ma_rad)) * m_zip_fraction;
            }

            // zippering is disabled
            else
            {
                m_z0calc = 0;
            }

            // catastrophe is enabled
            if (m_catastrophes_enabled)
            {
                // proportional catastrophe is enabled
                if (m_proportional_catastrophes)
                {
                    // magic angle > catastrophe start angle
                    if (m_magic_angle > m_cat_start_angle)
                    {
                        m_c0calc = (4.0 / PI) * 2.0 * m_induced_catastrophe_fraction
                                 * (1. - sin(ma_rad) + (ma_rad - cs_rad) * cos(ma_rad)) / (PI - 2. * cs_rad);
                    }

                    // magic angle <= catastrophe start angle
                    else
                    {
                        m_c0calc
                        = (4.0 / PI) * 2.0 * m_induced_catastrophe_fraction * (1.0 - sin(cs_rad)) / (PI - 2.0 * cs_rad);
                    }
                }

                // proportional catastrophe is disabled
                else
                {
                    m_c0calc = (4.0 / PI) * m_induced_catastrophe_fraction * cos(max(ma_rad, cs_rad));
                }
            }

            // catastrophe is disabled
            else
            {
                m_c0calc = 0;
            }

            // calculate x0, using the values of c0 and/or z0
            m_x0calc = (4.0 / PI) - m_z0calc - m_c0calc;
            break;

        // catastrophe first is requested
        case INT_CAT_FIRST:

            // catastrophe is enabled
            if (m_catastrophes_enabled)
            {
                // proportional catastrophe is enabled
                if (m_proportional_catastrophes)
                {
                    double temp = max(0., cs_rad);
                    m_c0calc = (4.0 / PI) * 2.0 * m_induced_catastrophe_fraction
                             * (1.0 - sin(temp) + (temp - cs_rad) * cos(temp)) / (PI - 2.0 * cs_rad);
                }

                // proportional catastrophe is disabled
                else
                {
                    m_c0calc = (4.0 / PI) * m_induced_catastrophe_fraction * cos(max(cs_rad, 0.));
                }
            }

            // zippering is enabled
            if (m_zippering_enabled)
            {
                m_z0calc = (4.0 / PI) * (1.0 - cos(ma_rad));

                // catastrophe is enabled and  (magic angle > catastrophe start angle)
                if (m_catastrophes_enabled && (cs_rad < ma_rad))
                {
                    // proportional catastrophe is enabled
                    if (m_proportional_catastrophes)
                    {
                        double left_boundary = max(0., cs_rad);
                        m_z0calc -= (4.0 / PI) * 2.0 * m_induced_catastrophe_fraction
                                  * ((left_boundary - cs_rad) * cos(left_boundary) - (ma_rad - cs_rad) * cos(ma_rad)
                                     - sin(left_boundary) + sin(ma_rad))
                                  / (PI - 2.0 * cs_rad);
                    }

                    // proportional catastrophe is disabled
                    else
                    {
                        m_z0calc -= (4.0 / PI) * m_induced_catastrophe_fraction * (cos(cs_rad) - cos(ma_rad));
                    }
                }
            }

            // zippering is disabled
            else
            {
                m_z0calc = 0;
            }

            // calculate x0, using the values of c0 and/or z0
            m_x0calc = (4.0 / PI) - m_z0calc - m_c0calc;
            break;

        default:
            cerr << "no function to calculate z0, c0, x0 etc. (interactionType = " << m_interaction_type << ")\n";
            return false;
    }

    // theoretial value of G
    double l0 = pow((2.0 * (-m_v_min + m_v_tm) * (m_v_plus - m_v_tm) * (m_v_plus - m_v_tm)) / (m_k_nuc * m_v_plus * (-m_v_min + m_v_plus)), 1. / 3.);
    double lavg_inv = m_k_res / (-m_v_min + m_v_tm) - m_k_cat / (m_v_plus - m_v_tm);
    double g_0 = l0 * lavg_inv;

    if (m_show_output == 1)
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
