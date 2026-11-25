#include "mesh.hpp"
#include "geometry.hpp"
#include "trajectory.hpp"
#include "region.hpp"
#include "system.hpp"
#include "mt_tip.hpp"
#include "parameters.hpp"
#include "linalg.hpp"

void Geometry::call_translator(SurfaceVector& s_vec, int region_index, int region_index_neigh)
{
    int index1(regions[region_index]->side_rev_map[region_index_neigh]);

    // use affine transformation of triangle to map a point from one region to other region
    Eigen::Vector2d v1(s_vec.x, s_vec.y), v2;
    v2 = regions[region_index]->side[index1].a * v1 + regions[region_index]->side[index1].b;
    s_vec.x = v2(0, 0);
    s_vec.y = v2(1, 0);

    s_vec.angle += regions[region_index]->side[index1].xy_rotation_angle;

    return;
}

TrajectoryVector Geometry::create_trajectory(const SurfaceVector& s_vec)
{
    // create a new trajectory with give trajectory vector
    return s_vec.region->insert_trajectory(s_vec);
}

TrajectoryVector Geometry::create_and_link_trajectory(
const SurfaceVector& new_base, Trajectory* oldtr, Direction olddir, double cos_angle, double p_cat)
{
    // first create a new trajectory vector (trajectory)
    TrajectoryVector t_vec = create_trajectory(new_base);

    // declare a dummy trajectory vector
    TrajectoryVector this_vec;

    // copy old trajectory vector
    this_vec.trajectory = oldtr;

    // old trajectory vector has direction forward, set the new trajectory vector as next
    if (olddir == ::FORWARD)
    {
        oldtr->next_tr = t_vec;
        oldtr->next_tr_cos_angle = cos_angle;
        oldtr->next_trp_cat = p_cat;

        // get the direction and position of the dummy trajectory vector
        this_vec.dir = BACKWARD;
        this_vec.pos = oldtr->length;
    }

    // old trajectory vector has direction backward, set the new trajectory vector as previous
    else
    {
        oldtr->prev_tr = t_vec;
        oldtr->prev_tr_cos_angle = cos_angle;
        oldtr->prev_trp_cat = p_cat;

        // get the direction and position of the dummy trajectory vector
        this_vec.dir = ::FORWARD;
        this_vec.pos = 0;
    }

    // new trajectory vector has direction forward, set the dummy trajectory vector as previous
    if (t_vec.dir == ::FORWARD)
    {
        t_vec.trajectory->prev_tr = this_vec;
        t_vec.trajectory->prev_tr_cos_angle = cos_angle;
        t_vec.trajectory->prev_trp_cat = p_cat;
    }

    // new trajectory vector has direction backward, set the dummy trajectory vector as next
    else
    {
        t_vec.trajectory->next_tr = this_vec;
        t_vec.trajectory->next_tr_cos_angle = cos_angle;
        t_vec.trajectory->next_trp_cat = p_cat;
    }

    return t_vec;
}

SurfaceVector Geometry::random_surface_vector()
{
    int region_selector{ 0 };
    double r_area{ 0.0 };
    int domain_type{ 0 };

    if (regions.size() > 1)
    {
        // nucleation on PPB
        if (system->p.pp_bk_nuc_fraction > 0)
        {
            double rand_gen = system->random_gen.rand();
            domain_type = (rand_gen < system->p.pp_bk_nuc_fraction) ? 1 : 2;
        }

        // nucleation outside PPB
        else
        {
            domain_type = 0;
        }

        // select a region
        r_area = system->random_gen.rand_dbl_exc(patch_area[domain_type]);

        if (r_area < 0.5 * patch_area[domain_type])
        {
            region_selector = 0;
            while (r_area > regions[regions_index[domain_type][region_selector]]->area)
            {
                r_area -= regions[regions_index[domain_type][region_selector++]]->area;

#ifdef DBG_ACID_TEST
                if (regionSelector == regions.size())
                {
                    cerr << "DBG/ASSERT: ERROR: Ran out of surfaces for random surface vector creation.\n";
                    exit(-1);
                }
#endif
            }
        }

        else
        {
            region_selector = regions_index[domain_type].size() - 1;
            r_area = patch_area[domain_type] - r_area;
            while (r_area > regions[regions_index[domain_type][region_selector]]->area)
            {
                r_area -= regions[regions_index[domain_type][region_selector--]]->area;

#ifdef DBG_ACID_TEST
                if (regionSelector == -1)
                {
                    cerr << "DBG/ASSERT: ERROR: Ran out of surfaces for random surface vector creation.\n";
                    exit(-1);
                }
#endif
            }
        }
    }

    return regions[regions_index[domain_type][region_selector]]->random_surface_vector();
}

bool Geometry::integrity_check()
{
    bool valid(true);
    int plus_grow_count(0);
    int plus_shrink_count(0);
    int minus_count(0);
    Trajectory* tptr(nullptr);

    for (size_t ridx = 0; ridx < regions.size(); ridx++)
    {
        plus_grow_count += regions[ridx]->growing_plus_tip_list.size();
        plus_shrink_count += regions[ridx]->shrinking_plus_tip_list.size();
        minus_count += regions[ridx]->minus_tip_list.size();
        tptr = regions[ridx]->trajectories.first();
        while (tptr != nullptr)
        {
            if (!tptr->integrity_check())
            {
                valid = false;
            }

            tptr = tptr->next();
        }
    }

    if (minus_count != plus_grow_count + plus_shrink_count)
    {
        cerr << "Integrity check failed: unequal number of plus and minus ends registered on regions\n";
        valid = false;
    }

    if (plus_grow_count != system->growing_mts.size())
    {
        cerr << "Integrity check failed: number of registered growing tips does not equal the number of growing "
                "microtubules.\n";
        valid = false;
    }

    if (plus_shrink_count != system->shrinking_mts.size())
    {
        cerr << "Integrity check failed: number of registered shrinking tips does not equal the number of shrinking "
                "microtubules.\n";
        valid = false;
    }

    int total_plus_grow_count(0);
    for (int d_tag = 0; d_tag <= system->p.face_number; d_tag++)
    {
        total_plus_grow_count += system->growing_tips_reg[d_tag];
    }

    if (plus_grow_count != total_plus_grow_count)
    {
        cerr << "Integrity check failed: number of growing tips does not equal cached values\n";
        valid = false;
    }

    return valid;
}

double Geometry::optical_length()
{
    double sum(0.0);

    // calculate total optical length present in the geometry
    for (size_t ridx = 0; ridx < regions.size(); ridx++)
    {
        sum += regions[ridx]->optical_length();
    }
    return sum;
}

int Geometry::trajectory_count()
{
    int sum(0);

    // calculate total number of trajectory present in the geometry
    for (size_t ridx = 0; ridx < regions.size(); ridx++)
    {
        sum += regions[ridx]->trajectories.size();
    }

    return sum;
}

TriMeshGeometry::TriMeshGeometry(System* s):
    Geometry(s, 0.0)
{
    // loading a triangle mesh
    pickup_shape(this);

    return;
}

void TriMeshGeometry::get_order_parameters(OrderParameters& op)
{
    OrderParametersRaw op_raw;

    // calculate order parameter components
    for (int eno = 0; eno < element_max; eno++)
    {
        static_cast<Cartesian*>(regions[eno])
        ->get_order_parameters_raw_flat(op_raw, regions[eno]->orientation, regions[eno]->area);

        // store local scalar order parameter
        op.local_order.push_back(op_raw.local_order);
        op.sv.push_back(op_raw.sv);
    }

    // normalizing the order parameter components
    if (op_raw.local_l > ZERO_CUTOFF)
    {
        op_raw.qxx /= op_raw.local_l;
        op_raw.qxy /= op_raw.local_l;
        op_raw.qxz /= op_raw.local_l;
        op_raw.qyy /= op_raw.local_l;
        op_raw.qyz /= op_raw.local_l;
        op_raw.qzz /= op_raw.local_l;
    }

    /// calculate Q(2) tensor order parameter
    op.r = op_raw.extract_r(op.rdirector, system->p.geometry);

    // normalize the order parameter director
    Eigen::Vector3d rvec(op.rdirector[0], op.rdirector[1], op.rdirector[2]);
    rvec /= rvec.norm();

    op.rdirector[0] = rvec(0, 0);
    op.rdirector[1] = rvec(1, 0);
    op.rdirector[2] = rvec(2, 0);

    op.c = fabs(op.rdirector[2])
             - 0.5 * sqrt(op.rdirector[0] * op.rdirector[0] + op.rdirector[1] * op.rdirector[1]);

    return;
}

void TriMeshGeometry::output_snapshot(ostream& out)
{
    // output the simulation results
    for (int eno = 0; eno < element_max; eno++)
    {
        regions[eno]->output_snapshot(out);
    }

    out << endl;

    return;
}

void TriMeshGeometry::output_order_heat_map(ostream& out, vector<double>& local_order, vector<Eigen::Vector3d>& sv)
{
    for (size_t eno = 0; eno < local_order.size(); eno++)
    {
        Eigen::Vector3d seg3_d;

        // Incircle radius of the triangle
        double sv_length = (2.0 * regions[eno]->area * regions[eno]->area) / regions[eno]->peri_meter;

        // get the region 3D rotation axis
        Eigen::Vector3d cross(regions[eno]->q.x(), regions[eno]->q.y(), regions[eno]->q.z());
        double cross_magnitude = cross.norm();

        Eigen::Quaternion<double> qr_1(0.0, 0.0, 0.0, 0.0), qr_2(0.0, 0.0, 0.0, 0.0);

        // region (and hence all the sgments in this region) need 3D rotation + translation
        if (cross_magnitude > 0.0)
        {
            // translate and rotate the start point of the segment
            qr_2.w() = 0;
            qr_2.x() = -1.0 * sv_length * sv[eno](0, 0) + regions[eno]->mid_point(0, 0);
            qr_2.y() = -1.0 * sv_length * sv[eno](1, 0) + regions[eno]->mid_point(1, 0);
            qr_2.z() = -1.0 * sv_length * sv[eno](2, 0) + regions[eno]->z_offset;

            qr_1 = regions[eno]->q.inverse() * qr_2 * regions[eno]->q;
            seg3_d << qr_1.x(), qr_1.y(), qr_1.z();
            seg3_d += object_cm;

            out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << ",";

            // translate and rotate the end point of the segment
            qr_2.w() = 0;
            qr_2.x() = sv_length * sv[eno](0, 0) + regions[eno]->mid_point(0, 0);
            qr_2.y() = sv_length * sv[eno](1, 0) + regions[eno]->mid_point(1, 0);
            qr_2.z() = sv_length * sv[eno](2, 0) + regions[eno]->z_offset;

            qr_1 = regions[eno]->q.inverse() * qr_2 * regions[eno]->q;
            seg3_d << qr_1.x(), qr_1.y(), qr_1.z();
            seg3_d += object_cm;

            out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << "," << regions[eno]->region_id
                << "@";
        }

        // need no rotation of the region (and hence all the sgments in this region)
        else
        {
            // translate the start point of the segment
            seg3_d << -1.0 * sv_length * sv[eno](0, 0) + regions[eno]->mid_point(0, 0),
            -1.0 * sv_length * sv[eno](1, 0) + regions[eno]->mid_point(1, 0),
            -1.0 * sv_length * sv[eno](2, 0) + regions[eno]->z_offset;
            seg3_d += object_cm;
            out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << ",";

            // translate the end point of the segment
            seg3_d << sv_length * sv[eno](0, 0) + regions[eno]->mid_point(0, 0),
            sv_length * sv[eno](1, 0) + regions[eno]->mid_point(1, 0),
            -1.0 * sv_length * sv[eno](2, 0) + regions[eno]->z_offset;
            seg3_d += object_cm;
            out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << "," << regions[eno]->region_id
                << "@";
        }
    }

    out << endl;

    return;
}

TrajectoryVector TriMeshGeometry::extend_trajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Triangle::extendTrajectory() called\n";
#endif

    // trajectory-director
    Eigen::Vector2d t_director(0.0, 0.0);

    SurfaceVector new_base(oldtr->base);

    // get copy of old trajectory region index
    int region_index(oldtr->base.region->region_id), region_index_neigh(0);

    // forward direction: move to the end of the trajectory
    if (dir == ::FORWARD)
    {
        new_base.x = oldtr->end_point[1].x;
        new_base.y = oldtr->end_point[1].y;

        // unit vector along the direction of old trajectory
        double t_dx = oldtr->end_point[1].x - oldtr->end_point[0].x;
        double t_dy = oldtr->end_point[1].y - oldtr->end_point[0].y;
        t_director << t_dx / oldtr->length, t_dy / oldtr->length;

        // get the neighbouring region index
        region_index_neigh = oldtr->end_point[1].next_element;
    }

    // backward direction: stay at the start of the trajectory
    if (dir == BACKWARD)
    {
        // unit vector along the direction of old trajectory
        double t_dx = oldtr->end_point[0].x - oldtr->end_point[1].x;
        double t_dy = oldtr->end_point[0].y - oldtr->end_point[1].y;
        t_director << t_dx / oldtr->length, t_dy / oldtr->length;

        // get the neighbouring region index
        region_index_neigh = oldtr->end_point[0].next_element;

        if (new_base.angle >= PI)
        {
            new_base.angle -= PI;
        }
        else
        {
            new_base.angle += PI;
        }
    }

    // get base coordinates for a fresh trajectory (to be created)
    call_translator(new_base, region_index, region_index_neigh);

    // get the new region for the fresh trajectory (to be created)
    new_base.region = regions[region_index_neigh];

    double cos_edge_angle(1.0), p_cat_reg(0.0);
    int s(oldtr->base.region->side_rev_map[region_index_neigh]);

    // cos-angle: between (encountered) triangle-edge and the trajectory-director
    double cos_bending_angle(t_director.dot(oldtr->base.region->side[s].dir2_d));

    // assign square of the cos-angle for this (encountered) triangle edge
    cos_edge_angle = cos_bending_angle * cos_bending_angle;

    // assign pCatReg for this for this (encountered) triangle edge
    p_cat_reg = oldtr->base.region->side[s].p_cat;

    // final correction to newBase angle -> [0 ... 2*PI]
    if (new_base.angle >= 2 * PI)
    {
        new_base.angle -= 2 * PI;
    }
    else if (new_base.angle < 0)
    {
        new_base.angle += 2 * PI;
    }

#ifdef DBG_GEOMETRY
    cout << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
#endif

    // create a fresh trajectory
    return create_and_link_trajectory(new_base, oldtr, dir, cos_edge_angle, p_cat_reg);
}

Cartesian::Cartesian(Geometry* g, double a):
    Region(g, a)
{
    // necessary to use many of its member functions
    return;
}

double Cartesian::intersection_angle(Trajectory* t1, Trajectory* t2)
{
#ifdef DBG_ACID_TEST
    if ((t2->base.angle < 0) || (t2->base.angle > 2 * PI))
    {
        cerr << "ERROR: angle out of bounds: " << t2->base.angle << ", region=" << RegionTypeText[t2->base.region->type]
             << "\n";
    }
    if ((t1->base.angle < 0) || (t1->base.angle > 2 * PI))
    {
        cerr << "ERROR: angle out of bounds: " << t1->base.angle << ", region=" << RegionTypeText[t1->base.region->type]
             << "\n";
    }
#endif

    // calculate collision angle between two interacting trajectories
    double angle(t2->base.angle - t1->base.angle);

    // angle is now in the interval -2PI..2PI
    if (angle < 0)
    {
        angle += 2 * PI; // to 0..2PI
    }

    if (angle > PI)
    {
        angle = 2 * PI - angle;
    }

    return (angle);
}

void Cartesian::get_order_parameters_raw_flat(OrderParametersRaw& op_r,
                                              vector<Eigen::Vector3d>& orientation,
                                              double area)
{
    double sin2(0.0), cos2(0.0);
    Eigen::Vector3d sv(0.0, 0.0, 0.0);
    double angle(0.0), length(0.0), local_length(0.);
    double qxx(0.0), qxy(0.0), qxz(0.0), qyy(0.0), qyz(0.0), qzz(0.0);
    double u1(0.0), u2(0.0), u3(0.0), order_local(0.0);
    Eigen::MatrixXd qf(3, 3), qc(3, 3), e(3, 3), o_l(2, 2);

    // get the rotation matrix
    e << orientation[0][0], orientation[1][0], orientation[2][0], orientation[0][1], orientation[1][1],
    orientation[2][1], orientation[0][2], orientation[1][2], orientation[2][2];

    Trajectory* tr(trajectories.first());

    while (tr != nullptr)
    {
        length = tr->segment_length();

        // add up all the segmnt lengths
        local_length += length;

        angle = tr->base.angle;

        // scalar (local) order parameter components
        sin2 += length * sin(2 * angle);
        cos2 += length * cos(2 * angle);

        // director for segment in local 2d-region
        u1 = cos(angle);
        u2 = sin(angle);
        u3 = 0.0;

        // Q(2) tensor (local) order parameter components
        qxx += length * (2 * u1 * u1);
        qxy += length * (2 * u1 * u2);
        qxz += length * (2 * u1 * u3);
        qyy += length * (2 * u2 * u2);
        qyz += length * (2 * u2 * u3);
        qzz += length * (2 * u3 * u3);

        tr = tr->next();
    }

    // only if trace of MT segment found
    if (local_length > 0)
    {
        qxx /= local_length;
        qxy /= local_length;
        qxz /= local_length;
        qyy /= local_length;
        qyz /= local_length;
        qzz /= local_length;

        sin2 /= local_length;
        cos2 /= local_length;

        // local order
        Eigen::Matrix3d matrix{ { qxx, qxy, qxz }, { qxy, qyy, qyz }, { qxz, qyz, qzz } };
        Eigen::Matrix3d evec_mat{ Eigen::Matrix3d::Zero() };
        Eigen::Vector3d e_val{ Eigen::Vector3d::Zero() };

        eigen_decomposition(matrix, evec_mat, e_val);

        int max_pos(0);

        if (e_val[1] > e_val[0])
        {
            max_pos = 1;
        }
        if (e_val[2] > e_val[max_pos])
        {
            max_pos = 2;
        }

        for (size_t i = 0; i < 3; ++i)
        {
            sv[i] = evec_mat(i, max_pos);
        }

        sv /= sv.norm();
        order_local = sqrt(sin2 * sin2 + cos2 * cos2) / area;

        // global order
        qxx -= 1;
        qyy -= 1;
        qf << qxx, qxy, qxz, qxy, qyy, qyz, qxz, qyz, qzz;

        // make tensor transformation
        qc = e * qf * e.transpose();

        // accumulate all the Q(2) tensors components from different regions
        op_r.qxx += local_length * qc(0, 0);
        op_r.qxy += local_length * qc(0, 1);
        op_r.qxz += local_length * qc(0, 2);
        op_r.qyy += local_length * qc(1, 1);
        op_r.qyz += local_length * qc(1, 2);
        op_r.qzz += local_length * qc(2, 2);

        // add up length of the segments from all the regions
        op_r.local_l += local_length;
    }

    op_r.local_order = order_local;
    op_r.sv = sv;
    return;
}

void Cartesian::output_snapshot(ostream& out)
{
    Cartesian::output_snapshot_offset(out, 0.0, 0.0);
    return;
}

void Cartesian::output_order_heat_map(ostream& out, vector<double>& local_order, vector<Eigen::Vector3d>& sv)
{
    return;
}

void Cartesian::output_snapshot_offset(ostream& out, double x_offset, double y_offset)
{
    Eigen::Vector3d seg3_d;

    // get the region 3D rotation axis
    Eigen::Vector3d cross(q.x(), q.y(), q.z());
    double cross_magnitude = cross.norm();

    Eigen::Quaternion<double> qr_1(0.0, 0.0, 0.0, 0.0), qr_2(0.0, 0.0, 0.0, 0.0);

    Trajectory* tr;
    SurfaceVector svec;
    list<Segment*>::iterator seg;

    tr = trajectories.first();

    while (tr != nullptr)
    {
        svec = tr->base;

        // record segments in the snapshots
        seg = tr->segments.begin();
        while (seg != tr->segments.end())
        {
            svec = tr->base;

            // move to the starting end of the segment
            translate_vector(svec, (**seg).start);

            // for backward direction retrieve actual base angle (during nucleation)
            if ((**seg).dir == BACKWARD)
            {
                svec.angle -= PI;
                if (svec.angle < 0)
                {
                    svec.angle += 2 * PI;
                }
            }

            // region (and hence all the sgments in this region) need 3D rotation + translation
            if (cross_magnitude > 0.0)
            {
                // translate and rotate the start point of the segment
                qr_2.w() = 0;
                qr_2.x() = svec.x + mid_point(0, 0);
                qr_2.y() = svec.y + mid_point(1, 0);
                qr_2.z() = z_offset;

                qr_1 = q.inverse() * qr_2 * q;
                seg3_d << qr_1.x(), qr_1.y(), qr_1.z();
                seg3_d += geometry->object_cm;

                out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << ",";

                // translate and rotate the end point of the segment
                qr_2.w() = 0;
                qr_2.x() = svec.x + mid_point(0, 0) + (**seg).length() * cos(svec.angle);
                qr_2.y() = svec.y + mid_point(1, 0) + (**seg).length() * sin(svec.angle);
                qr_2.z() = z_offset;

                qr_1 = q.inverse() * qr_2 * q;
                seg3_d << qr_1.x(), qr_1.y(), qr_1.z();
                seg3_d += geometry->object_cm;

                out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << "," << region_id << "@";
            }

            // need no rotation of the region (and hence all the sgments in this region)
            else
            {
                // translate the start point of the segment
                seg3_d << svec.x + mid_point(0, 0), svec.y + mid_point(1, 0), z_offset;
                seg3_d += geometry->object_cm;
                out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << ",";

                // translate the end point of the segment
                seg3_d << svec.x + mid_point(0, 0) + (**seg).length() * cos(svec.angle),
                svec.y + mid_point(1, 0) + (**seg).length() * sin(svec.angle), z_offset;
                seg3_d += geometry->object_cm;
                out << seg3_d(0, 0) << "," << seg3_d(1, 0) << "," << seg3_d(2, 0) << "," << region_id << "@";
            }
            seg++;
        }

        tr = tr->next();
    }

    return;
}

void Cartesian::translate_vector(SurfaceVector& s_vec, const double dist)
{
    // this can be used to move a point along a trajectory (used in making snap-shots)
    s_vec.x += cos(s_vec.angle) * dist;
    s_vec.y += sin(s_vec.angle) * dist;
    return;
}

void Cartesian::make_intersection_list(Trajectory* tr1)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Cartesian::make_intersection_list() called\n";
#endif

    double cos_diff(0.0);
    double denominator(0.0);
    double temp1(0.0), temp2(0.0);
    double bp1(0.0), bp2(0.0);
    double cos1(0.0), cos2(0.0), sin1(0.0), sin2(0.0);

    // take cos/sin of first trajectory base angle
    cos1 = cos(tr1->base.angle);
    sin1 = sin(tr1->base.angle);

    Intersection is_this;
    Intersection is_that;
    IntersectionItr this_itr;
    IntersectionItr that_itr;

    is_this.occupancy = 0;
    is_that.occupancy = 0;
    is_that.other_trajectory = tr1;

    Trajectory* tr2(trajectories.first());

    while (tr2 != nullptr)
    {
        // imagine a triangle with vertices sitting at: (a) one at the trajectory intersection point and (b) other two
        // at the base of the two trajectories
        if (tr1 != tr2)
        {
            // take cos/sin of first trajectory base angle
            cos2 = cos(tr2->base.angle);
            sin2 = sin(tr2->base.angle);

            cos_diff = cos1 * cos2 + sin1 * sin2;
            denominator = 1 - cos_diff * cos_diff;

            // find the intersection point distance on the trajectory
            if (denominator > ZERO_CUTOFF)
            {
                denominator = 1 / denominator;
                temp1 = (tr2->base.x - tr1->base.x) * denominator;
                temp2 = (tr2->base.y - tr1->base.y) * denominator;

                bp1 = temp1 * (cos1 - cos_diff * cos2) + temp2 * (sin1 - cos_diff * sin2);

                // intersetion located within the range of trajectory length
                if (!((bp1 < 0) || (bp1 > tr1->length)))
                {
                    bp2 = -temp1 * (cos2 - cos_diff * cos1) - temp2 * (sin2 - cos_diff * sin1);
#ifdef DBG_ASSERT
                    if ((bp2 < 0) || (bp2 > tr2->length))
                    {
                        cout << "ERROR: DBG/extra check: incompatible results for line intersection.\n";
                        cout << "tr1: " << bp1 << " out of " << tr1->length << ", base=" << tr1->base.x << ","
                             << tr1->base.y << " angle=" << tr1->base.angle << "\n";
                        cout << "tr2: " << bp2 << " out of " << tr2->length << ", base=" << tr2->base.x << ","
                             << tr2->base.y << " angle=" << tr2->base.angle << "\n";
                        exit(-2);
                    }
#endif

                    // first insert 'isThat': occupancy and other are still incorrect
                    that_itr = tr2->intersections.insert(pair<double, Intersection>(bp2, is_that));
                    is_this.mirror = that_itr;
                    is_this.other_trajectory = tr2;
                    this_itr = tr1->intersections.insert(pair<double, Intersection>(bp1, is_this));

                    // update the reference in the other half of the intersection
                    that_itr->second.mirror = this_itr;

                    // call the update routine of the other trajectory, this will also correct the occupancy number
                    tr2->new_intersection(that_itr);
                }
            }
        }

        tr2 = tr2->next();
    }
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: " << tr1->intersections.size() - 1 << " intersections\n";
#endif

    return;
}

Triangle::Triangle(double element_area, Geometry* g):
    Cartesian(g, element_area)
{
    // necessary to create a full geometry, made with many triangles
    return;
}

SurfaceVector Triangle::random_surface_vector()
{
    SurfaceVector v;
    double r1(0.0), r2(0.0);

    // get system pointer
    System* s(geometry->system);

    // random number between (0,1)
    r1 = s->random_gen.rand_dbl_exc(1.0);
    r2 = s->random_gen.rand_dbl_exc(1.0);

    // store sqrt(r1) to use it many times, hence avoid expensive call of sqrt()
    double sqr1 = sqrt(r1);

    // generate uniformly and isotropically distributed nucleation points on each triangle
    v.x
    = (1 - sqr1) * vertices[0](0, 0) + (sqr1 * (1 - r2)) * vertices[1](0, 0) + (sqr1 * r2) * vertices[2](0, 0);
    v.y
    = (1 - sqr1) * vertices[0](1, 0) + (sqr1 * (1 - r2)) * vertices[1](1, 0) + (sqr1 * r2) * vertices[2](1, 0);

    // make the unifor distribution also isotropic
    v.angle = s->random_gen.rand_exc(2 * PI);

    v.region = this;

    return v;
}

void Triangle::get_trajectory_coordinates(SurfaceVector& s_vec,
                                          double& total_length,
                                          vector<PointATedge>& end_point,
                                          TrajectoryVector& t_vec)
{
    // copy old base coordinate
    double old_x(s_vec.x);
    double old_y(s_vec.y);

    t_vec.trajectory = nullptr;

    double acos(cos(s_vec.angle));
    double asin(sin(s_vec.angle));

    // direction and base angle of the trajectory
    if (abs(acos) < ZERO_CUTOFF)
    {
        if (s_vec.angle < PI)
        {
            t_vec.dir = ::FORWARD;
        }
        else
        {
            t_vec.dir = BACKWARD;
            s_vec.angle -= PI;
        }
    }

    else if (acos < 0)
    {
        s_vec.angle += PI;
        t_vec.dir = BACKWARD;
        if (s_vec.angle > 2 * PI)
        {
            s_vec.angle -= 2 * PI;
        }
    }

    else
    {
        t_vec.dir = ::FORWARD;
    }

    // use the triangle perimeter to create a trajectory, this will to make sure to have exactly two intersections
    Eigen::Vector2d v1(s_vec.x + peri_meter * acos, s_vec.y + peri_meter * asin);
    Eigen::Vector2d v2(s_vec.x - peri_meter * acos, s_vec.y - peri_meter * asin);

    // intersection points of the trajectory with the triangle edges
    int encounter(0), neno(0);
    for (int i = 0; i < 3; i++)
    {
        if (encounter == 2)
        {
            break;
        }

        Eigen::Vector2d u1;
        u1 = vertices[side[i].dir[0]];

        Eigen::Vector2d u2;
        u2 = vertices[side[i].dir[1]];

        Eigen::Vector2d u;
        u = u2 - u1;

        Eigen::Vector2d v;
        v = v2 - v1;

        Eigen::Vector2d w;
        w = u1 - v1;

        double sval(0.0);
        sval = (v(1, 0) * w(0, 0) - v(0, 0) * w(1, 0)) / (v(0, 0) * u(1, 0) - v(1, 0) * u(0, 0));

        if ((sval >= 0.0) && (sval <= 1.0))
        {
            encounter += 1;
            neno = side_map[i];

            // get triangle-edge to trajectory intersection point
            end_point.push_back(
            PointATedge(u1(0, 0) + sval * (u2(0, 0) - u1(0, 0)), u1(1, 0) + sval * (u2(1, 0) - u1(1, 0)), 0.0, neno));
        }
    }

    // base coordinate (x_base < x_other) of the trajectory
    if (end_point[0].x > end_point[1].x)
    {
        reverse(end_point.begin(), end_point.end());
    }

    s_vec.x = end_point[0].x;
    s_vec.y = end_point[0].y;

    // nucleation position (distance between old base and nucleation point) on the trajectory
    double dx(s_vec.x - old_x), dy(s_vec.y - old_y);
    t_vec.pos = sqrt(dx * dx + dy * dy);

    // length of the trajectory (distance between the end points of the trajectory)
    dx = end_point[1].x - end_point[0].x;
    dy = end_point[1].y - end_point[0].y;
    total_length = sqrt(dx * dx + dy * dy);

    // eacch trajectory should have exacty two intersection point with exactly two triangle edges
    if (encounter != 2)
    {
        cout << "Problem with the number of trajectory-triangle intersections (should be == 2) ---> " << encounter
             << endl;
        exit(-1);
    }

    return;
}
