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
    int Index1(regions[region_index]->sideRevMap[region_index_neigh]);

    // use affine transformation of triangle to map a point from one region to other region
    Eigen::Vector2d V1(s_vec.x, s_vec.y), V2;
    V2 = regions[region_index]->side[Index1].A * V1 + regions[region_index]->side[Index1].b;
    s_vec.x = V2(0, 0);
    s_vec.y = V2(1, 0);

    s_vec.angle += regions[region_index]->side[Index1].xyRotationAngle;

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
    TrajectoryVector tVec = create_trajectory(new_base);

    // declare a dummy trajectory vector
    TrajectoryVector thisVec;

    // copy old trajectory vector
    thisVec.trajectory = oldtr;

    // old trajectory vector has direction forward, set the new trajectory vector as next
    if (olddir == ::FORWARD)
    {
        oldtr->nextTr = tVec;
        oldtr->nextTrCosAngle = cos_angle;
        oldtr->nextTrpCat = p_cat;

        // get the direction and position of the dummy trajectory vector
        thisVec.dir = BACKWARD;
        thisVec.pos = oldtr->length;
    }

    // old trajectory vector has direction backward, set the new trajectory vector as previous
    else
    {
        oldtr->prevTr = tVec;
        oldtr->prevTrCosAngle = cos_angle;
        oldtr->prevTrpCat = p_cat;

        // get the direction and position of the dummy trajectory vector
        thisVec.dir = ::FORWARD;
        thisVec.pos = 0;
    }

    // new trajectory vector has direction forward, set the dummy trajectory vector as previous
    if (tVec.dir == ::FORWARD)
    {
        tVec.trajectory->prevTr = thisVec;
        tVec.trajectory->prevTrCosAngle = cos_angle;
        tVec.trajectory->prevTrpCat = p_cat;
    }

    // new trajectory vector has direction backward, set the dummy trajectory vector as next
    else
    {
        tVec.trajectory->nextTr = thisVec;
        tVec.trajectory->nextTrCosAngle = cos_angle;
        tVec.trajectory->nextTrpCat = p_cat;
    }

    return tVec;
}

SurfaceVector Geometry::random_surface_vector()
{
    int regionSelector{ 0 };
    double rArea{ 0.0 };
    int domainType{ 0 };

    if (regions.size() > 1)
    {
        // nucleation on PPB
        if (system->p.PPBkNucFraction > 0)
        {
            double randGen = system->randomGen.rand();
            domainType = (randGen < system->p.PPBkNucFraction) ? 1 : 2;
        }

        // nucleation outside PPB
        else
        {
            domainType = 0;
        }

        // select a region
        rArea = system->randomGen.rand_dbl_exc(patchArea[domainType]);

        if (rArea < 0.5 * patchArea[domainType])
        {
            regionSelector = 0;
            while (rArea > regions[RegionsIndex[domainType][regionSelector]]->area)
            {
                rArea -= regions[RegionsIndex[domainType][regionSelector++]]->area;

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
            regionSelector = RegionsIndex[domainType].size() - 1;
            rArea = patchArea[domainType] - rArea;
            while (rArea > regions[RegionsIndex[domainType][regionSelector]]->area)
            {
                rArea -= regions[RegionsIndex[domainType][regionSelector--]]->area;

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

    return regions[RegionsIndex[domainType][regionSelector]]->random_surface_vector();
}

bool Geometry::integrity_check()
{
    bool valid(true);
    int plusGrowCount(0);
    int plusShrinkCount(0);
    int minusCount(0);
    Trajectory* tptr(nullptr);

    for (size_t ridx = 0; ridx < regions.size(); ridx++)
    {
        plusGrowCount += regions[ridx]->growingPlusTipList.size();
        plusShrinkCount += regions[ridx]->shrinkingPlusTipList.size();
        minusCount += regions[ridx]->minusTipList.size();
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

    if (minusCount != plusGrowCount + plusShrinkCount)
    {
        cerr << "Integrity check failed: unequal number of plus and minus ends registered on regions\n";
        valid = false;
    }

    if (plusGrowCount != system->growing_mts.size())
    {
        cerr << "Integrity check failed: number of registered growing tips does not equal the number of growing "
                "microtubules.\n";
        valid = false;
    }

    if (plusShrinkCount != system->shrinking_mts.size())
    {
        cerr << "Integrity check failed: number of registered shrinking tips does not equal the number of shrinking "
                "microtubules.\n";
        valid = false;
    }

    int totalPlusGrowCount(0);
    for (int dTag = 0; dTag <= system->p.faceNumber; dTag++)
    {
        totalPlusGrowCount += system->growingTipsReg[dTag];
    }

    if (plusGrowCount != totalPlusGrowCount)
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
    OrderParametersRaw opRaw;

    // calculate order parameter components
    for (int eno = 0; eno < elementMax; eno++)
    {
        static_cast<Cartesian*>(regions[eno])
        ->get_order_parameters_raw_flat(opRaw, regions[eno]->orientation, regions[eno]->area);

        // store local scalar order parameter
        op.localOrder.push_back(opRaw.localOrder);
        op.Sv.push_back(opRaw.Sv);
    }

    // normalizing the order parameter components
    if (opRaw.localL > ZERO_CUTOFF)
    {
        opRaw.Qxx /= opRaw.localL;
        opRaw.Qxy /= opRaw.localL;
        opRaw.Qxz /= opRaw.localL;
        opRaw.Qyy /= opRaw.localL;
        opRaw.Qyz /= opRaw.localL;
        opRaw.Qzz /= opRaw.localL;
    }

    /// calculate Q(2) tensor order parameter
    op.R = opRaw.extract_r(op.Rdirector, system->p.geometry);

    // normalize the order parameter director
    Eigen::Vector3d Rvec(op.Rdirector[0], op.Rdirector[1], op.Rdirector[2]);
    Rvec /= Rvec.norm();

    op.Rdirector[0] = Rvec(0, 0);
    op.Rdirector[1] = Rvec(1, 0);
    op.Rdirector[2] = Rvec(2, 0);

    op.C = fabs(op.Rdirector[2]) - 0.5 * sqrt(op.Rdirector[0] * op.Rdirector[0] + op.Rdirector[1] * op.Rdirector[1]);

    return;
}

void TriMeshGeometry::output_snapshot(ostream& out)
{
    // output the simulation results
    for (int eno = 0; eno < elementMax; eno++)
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
        Eigen::Vector3d seg3D;

        // Incircle radius of the triangle
        double svLength = (2.0 * regions[eno]->area * regions[eno]->area) / regions[eno]->periMeter;

        // get the region 3D rotation axis
        Eigen::Vector3d cross(regions[eno]->Q.x(), regions[eno]->Q.y(), regions[eno]->Q.z());
        double crossMagnitude = cross.norm();

        Eigen::Quaternion<double> qr_1(0.0, 0.0, 0.0, 0.0), qr_2(0.0, 0.0, 0.0, 0.0);

        // region (and hence all the sgments in this region) need 3D rotation + translation
        if (crossMagnitude > 0.0)
        {
            // translate and rotate the start point of the segment
            qr_2.w() = 0;
            qr_2.x() = -1.0 * svLength * sv[eno](0, 0) + regions[eno]->midPoint(0, 0);
            qr_2.y() = -1.0 * svLength * sv[eno](1, 0) + regions[eno]->midPoint(1, 0);
            qr_2.z() = -1.0 * svLength * sv[eno](2, 0) + regions[eno]->zOffset;

            qr_1 = regions[eno]->Q.inverse() * qr_2 * regions[eno]->Q;
            seg3D << qr_1.x(), qr_1.y(), qr_1.z();
            seg3D += objectCM;

            out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << ",";

            // translate and rotate the end point of the segment
            qr_2.w() = 0;
            qr_2.x() = svLength * sv[eno](0, 0) + regions[eno]->midPoint(0, 0);
            qr_2.y() = svLength * sv[eno](1, 0) + regions[eno]->midPoint(1, 0);
            qr_2.z() = svLength * sv[eno](2, 0) + regions[eno]->zOffset;

            qr_1 = regions[eno]->Q.inverse() * qr_2 * regions[eno]->Q;
            seg3D << qr_1.x(), qr_1.y(), qr_1.z();
            seg3D += objectCM;

            out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << "," << regions[eno]->regionId << "@";
        }

        // need no rotation of the region (and hence all the sgments in this region)
        else
        {
            // translate the start point of the segment
            seg3D << -1.0 * svLength * sv[eno](0, 0) + regions[eno]->midPoint(0, 0),
            -1.0 * svLength * sv[eno](1, 0) + regions[eno]->midPoint(1, 0),
            -1.0 * svLength * sv[eno](2, 0) + regions[eno]->zOffset;
            seg3D += objectCM;
            out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << ",";

            // translate the end point of the segment
            seg3D << svLength * sv[eno](0, 0) + regions[eno]->midPoint(0, 0),
            svLength * sv[eno](1, 0) + regions[eno]->midPoint(1, 0),
            -1.0 * svLength * sv[eno](2, 0) + regions[eno]->zOffset;
            seg3D += objectCM;
            out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << "," << regions[eno]->regionId << "@";
        }
    }

    out << endl;

    return;
}

TrajectoryVector TriMeshGeometry::extend_trajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Triangle::extend_trajectory() called\n";
#endif

    // trajectory-director
    Eigen::Vector2d tDirector(0.0, 0.0);

    SurfaceVector newBase(oldtr->base);

    // get copy of old trajectory region index
    int regionIndex(oldtr->base.region->regionId), regionIndexNeigh(0);

    // forward direction: move to the end of the trajectory
    if (dir == ::FORWARD)
    {
        newBase.x = oldtr->endPoint[1].x;
        newBase.y = oldtr->endPoint[1].y;

        // unit vector along the direction of old trajectory
        double tDx = oldtr->endPoint[1].x - oldtr->endPoint[0].x;
        double tDy = oldtr->endPoint[1].y - oldtr->endPoint[0].y;
        tDirector << tDx / oldtr->length, tDy / oldtr->length;

        // get the neighbouring region index
        regionIndexNeigh = oldtr->endPoint[1].nextElement;
    }

    // backward direction: stay at the start of the trajectory
    if (dir == BACKWARD)
    {
        // unit vector along the direction of old trajectory
        double tDx = oldtr->endPoint[0].x - oldtr->endPoint[1].x;
        double tDy = oldtr->endPoint[0].y - oldtr->endPoint[1].y;
        tDirector << tDx / oldtr->length, tDy / oldtr->length;

        // get the neighbouring region index
        regionIndexNeigh = oldtr->endPoint[0].nextElement;

        if (newBase.angle >= PI)
        {
            newBase.angle -= PI;
        }
        else
        {
            newBase.angle += PI;
        }
    }

    // get base coordinates for a fresh trajectory (to be created)
    call_translator(newBase, regionIndex, regionIndexNeigh);

    // get the new region for the fresh trajectory (to be created)
    newBase.region = regions[regionIndexNeigh];

    double cosEdgeAngle(1.0), pCatReg(0.0);
    int s(oldtr->base.region->sideRevMap[regionIndexNeigh]);

    // cos-angle: between (encountered) triangle-edge and the trajectory-director
    double cosBendingAngle(tDirector.dot(oldtr->base.region->side[s].dir2D));

    // assign square of the cos-angle for this (encountered) triangle edge
    cosEdgeAngle = cosBendingAngle * cosBendingAngle;

    // assign pCatReg for this for this (encountered) triangle edge
    pCatReg = oldtr->base.region->side[s].pCat;

    // final correction to newBase angle -> [0 ... 2*PI]
    if (newBase.angle >= 2 * PI)
    {
        newBase.angle -= 2 * PI;
    }
    else if (newBase.angle < 0)
    {
        newBase.angle += 2 * PI;
    }

#ifdef DBG_GEOMETRY
    cout << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
#endif

    // create a fresh trajectory
    return create_and_link_trajectory(newBase, oldtr, dir, cosEdgeAngle, pCatReg);
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
    double angle(0.0), length(0.0), localLength(0.);
    double qxx(0.0), qxy(0.0), qxz(0.0), qyy(0.0), qyz(0.0), qzz(0.0);
    double u1(0.0), u2(0.0), u3(0.0), orderLocal(0.0);
    Eigen::MatrixXd QF(3, 3), QC(3, 3), e(3, 3), O_l(2, 2);

    // get the rotation matrix
    e << orientation[0][0], orientation[1][0], orientation[2][0], orientation[0][1], orientation[1][1],
    orientation[2][1], orientation[0][2], orientation[1][2], orientation[2][2];

    Trajectory* tr(trajectories.first());

    while (tr != nullptr)
    {
        length = tr->segment_length();

        // add up all the segmnt lengths
        localLength += length;

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
    if (localLength > 0)
    {
        qxx /= localLength;
        qxy /= localLength;
        qxz /= localLength;
        qyy /= localLength;
        qyz /= localLength;
        qzz /= localLength;

        sin2 /= localLength;
        cos2 /= localLength;

        // local order
        Eigen::Matrix3d matrix{ { qxx, qxy, qxz }, { qxy, qyy, qyz }, { qxz, qyz, qzz } };
        Eigen::Matrix3d evecMat{ Eigen::Matrix3d::Zero() };
        Eigen::Vector3d eVal{ Eigen::Vector3d::Zero() };

        eigen_decomposition(matrix, evecMat, eVal);

        int maxPos(0);

        if (eVal[1] > eVal[0])
        {
            maxPos = 1;
        }
        if (eVal[2] > eVal[maxPos])
        {
            maxPos = 2;
        }

        for (size_t i = 0; i < 3; ++i)
        {
            sv[i] = evecMat(i, maxPos);
        }

        sv /= sv.norm();
        orderLocal = sqrt(sin2 * sin2 + cos2 * cos2) / area;

        // global order
        qxx -= 1;
        qyy -= 1;
        QF << qxx, qxy, qxz, qxy, qyy, qyz, qxz, qyz, qzz;

        // make tensor transformation
        QC = e * QF * e.transpose();

        // accumulate all the Q(2) tensors components from different regions
        op_r.Qxx += localLength * QC(0, 0);
        op_r.Qxy += localLength * QC(0, 1);
        op_r.Qxz += localLength * QC(0, 2);
        op_r.Qyy += localLength * QC(1, 1);
        op_r.Qyz += localLength * QC(1, 2);
        op_r.Qzz += localLength * QC(2, 2);

        // add up length of the segments from all the regions
        op_r.localL += localLength;
    }

    op_r.localOrder = orderLocal;
    op_r.Sv = sv;
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
    Eigen::Vector3d seg3D;

    // get the region 3D rotation axis
    Eigen::Vector3d cross(Q.x(), Q.y(), Q.z());
    double crossMagnitude = cross.norm();

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
            if (crossMagnitude > 0.0)
            {
                // translate and rotate the start point of the segment
                qr_2.w() = 0;
                qr_2.x() = svec.x + midPoint(0, 0);
                qr_2.y() = svec.y + midPoint(1, 0);
                qr_2.z() = zOffset;

                qr_1 = Q.inverse() * qr_2 * Q;
                seg3D << qr_1.x(), qr_1.y(), qr_1.z();
                seg3D += geometry->objectCM;

                out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << ",";

                // translate and rotate the end point of the segment
                qr_2.w() = 0;
                qr_2.x() = svec.x + midPoint(0, 0) + (**seg).length() * cos(svec.angle);
                qr_2.y() = svec.y + midPoint(1, 0) + (**seg).length() * sin(svec.angle);
                qr_2.z() = zOffset;

                qr_1 = Q.inverse() * qr_2 * Q;
                seg3D << qr_1.x(), qr_1.y(), qr_1.z();
                seg3D += geometry->objectCM;

                out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << "," << regionId << "@";
            }

            // need no rotation of the region (and hence all the sgments in this region)
            else
            {
                // translate the start point of the segment
                seg3D << svec.x + midPoint(0, 0), svec.y + midPoint(1, 0), zOffset;
                seg3D += geometry->objectCM;
                out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << ",";

                // translate the end point of the segment
                seg3D << svec.x + midPoint(0, 0) + (**seg).length() * cos(svec.angle),
                svec.y + midPoint(1, 0) + (**seg).length() * sin(svec.angle), zOffset;
                seg3D += geometry->objectCM;
                out << seg3D(0, 0) << "," << seg3D(1, 0) << "," << seg3D(2, 0) << "," << regionId << "@";
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
    cout << "DBG/GEOMETRY: Cartesian::makeIntersectionList() called\n";
#endif

    double cosDiff(0.0);
    double denominator(0.0);
    double temp1(0.0), temp2(0.0);
    double bp1(0.0), bp2(0.0);
    double cos1(0.0), cos2(0.0), sin1(0.0), sin2(0.0);

    // take cos/sin of first trajectory base angle
    cos1 = cos(tr1->base.angle);
    sin1 = sin(tr1->base.angle);

    Intersection isThis;
    Intersection isThat;
    IntersectionItr thisItr;
    IntersectionItr thatItr;

    isThis.occupancy = 0;
    isThat.occupancy = 0;
    isThat.otherTrajectory = tr1;

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

            cosDiff = cos1 * cos2 + sin1 * sin2;
            denominator = 1 - cosDiff * cosDiff;

            // find the intersection point distance on the trajectory
            if (denominator > ZERO_CUTOFF)
            {
                denominator = 1 / denominator;
                temp1 = (tr2->base.x - tr1->base.x) * denominator;
                temp2 = (tr2->base.y - tr1->base.y) * denominator;

                bp1 = temp1 * (cos1 - cosDiff * cos2) + temp2 * (sin1 - cosDiff * sin2);

                // intersetion located within the range of trajectory length
                if (!((bp1 < 0) || (bp1 > tr1->length)))
                {
                    bp2 = -temp1 * (cos2 - cosDiff * cos1) - temp2 * (sin2 - cosDiff * sin1);
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
                    thatItr = tr2->intersections.insert(pair<double, Intersection>(bp2, isThat));
                    isThis.mirror = thatItr;
                    isThis.otherTrajectory = tr2;
                    thisItr = tr1->intersections.insert(pair<double, Intersection>(bp1, isThis));

                    // update the reference in the other half of the intersection
                    thatItr->second.mirror = thisItr;

                    // call the update routine of the other trajectory, this will also correct the occupancy number
                    tr2->new_intersection(thatItr);
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
    r1 = s->randomGen.rand_dbl_exc(1.0);
    r2 = s->randomGen.rand_dbl_exc(1.0);

    // store sqrt(r1) to use it many times, hence avoid expensive call of sqrt()
    double sqr1 = sqrt(r1);

    // generate uniformly and isotropically distributed nucleation points on each triangle
    v.x = (1 - sqr1) * vertices[0](0, 0) + (sqr1 * (1 - r2)) * vertices[1](0, 0) + (sqr1 * r2) * vertices[2](0, 0);
    v.y = (1 - sqr1) * vertices[0](1, 0) + (sqr1 * (1 - r2)) * vertices[1](1, 0) + (sqr1 * r2) * vertices[2](1, 0);

    // make the unifor distribution also isotropic
    v.angle = s->randomGen.rand_exc(2 * PI);

    v.region = this;

    return v;
}

void Triangle::get_trajectory_coordinates(SurfaceVector& s_vec,
                                          double& total_length,
                                          vector<PointATedge>& end_point,
                                          TrajectoryVector& t_vec)
{
    // copy old base coordinate
    double oldX(s_vec.x);
    double oldY(s_vec.y);

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
    Eigen::Vector2d v1(s_vec.x + periMeter * acos, s_vec.y + periMeter * asin);
    Eigen::Vector2d v2(s_vec.x - periMeter * acos, s_vec.y - periMeter * asin);

    // intersection points of the trajectory with the triangle edges
    int encounter(0), Neno(0);
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
            Neno = sideMap[i];

            // get triangle-edge to trajectory intersection point
            end_point.push_back(
            PointATedge(u1(0, 0) + sval * (u2(0, 0) - u1(0, 0)), u1(1, 0) + sval * (u2(1, 0) - u1(1, 0)), 0.0, Neno));
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
    double dx(s_vec.x - oldX), dy(s_vec.y - oldY);
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
