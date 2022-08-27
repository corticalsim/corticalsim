#include "corticalSimReal.h"
#include "meshImport.h"

TriMeshGeometry::TriMeshGeometry(System* s):Geometry(s,0.0)
{
	// loading a triangle mesh 
    	pickup_shape(this);

    	return;
}

void TriMeshGeometry::getOrderParameters(OrderParameters& op)
{
	OrderParametersRaw opRaw;

	// calculate order parameter components
    	for(int eno=0;eno<elementMax;eno++)
    	{
    		static_cast<Cartesian*>(regions[eno])->getOrderParametersRawFlat(opRaw,regions[eno]->orientation,regions[eno]->area);

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
    	op.R = opRaw.extractR(op.Rdirector,system->p.geometry);
    
    	// normalize the order parameter director
        Vector3d Rvec(op.Rdirector[0],op.Rdirector[1],op.Rdirector[2]);
	Rvec/= Rvec.norm();

    	op.Rdirector[0] = Rvec(0,0);
	op.Rdirector[1] = Rvec(1,0);
	op.Rdirector[2] = Rvec(2,0);

        op.C = fabs(op.Rdirector[2])-0.5*sqrt(op.Rdirector[0]*op.Rdirector[0]+op.Rdirector[1]*op.Rdirector[1]);

    return;
}

void TriMeshGeometry::outputSnapshot(ostream& out)
{
    // output the simulation results 
    for(int eno = 0;eno< elementMax; eno++)
    regions[eno]->outputSnapshot(out);
   
    out << endl;

    return;
}

void TriMeshGeometry::outputOrderHeatMap(ostream& out,vector<double>& localOrder,vector<Vector3d>& sv)
{
	for(int eno = 0; eno < localOrder.size(); eno++)
	{
        	Vector3d seg3D;

		// Incircle radius of the triangle 
		double svLength = (2.0*regions[eno]->area*regions[eno]->area)/regions[eno]->periMeter;

        	// get the region 3D rotation axis 
		Vector3d cross(regions[eno]->Q.x(),regions[eno]->Q.y(),regions[eno]->Q.z());
     		double crossMagnitude = cross.norm();  
  
		Quaternion<double> qr(0.0,0.0,0.0,0.0),QR(0.0,0.0,0.0,0.0);


                // region (and hence all the sgments in this region) need 3D rotation + translation
		if(crossMagnitude>0.0)
		{
                        // translate and rotate the start point of the segment 
			QR.w()= 0;
			QR.x() = -1.0*svLength*sv[eno](0,0)+regions[eno]->midPoint(0,0);
			QR.y() = -1.0*svLength*sv[eno](1,0)+regions[eno]->midPoint(1,0);
			QR.z() = -1.0*svLength*sv[eno](2,0)+regions[eno]->zOffset;

			qr = regions[eno]->Q.inverse() * QR * regions[eno]->Q;
			seg3D << qr.x(),qr.y(),qr.z();
			seg3D+=objectCM;

			out << seg3D(0,0) << "," <<seg3D(1,0)<< "," << seg3D(2,0) << ",";

                        // translate and rotate the end point of the segment
			QR.w()= 0;
			QR.x() = svLength*sv[eno](0,0)+regions[eno]->midPoint(0,0);
			QR.y() = svLength*sv[eno](1,0)+regions[eno]->midPoint(1,0);
			QR.z() = svLength*sv[eno](2,0)+regions[eno]->zOffset;

			qr = regions[eno]->Q.inverse() * QR * regions[eno]->Q;
			seg3D << qr.x(),qr.y(),qr.z();
			seg3D+=objectCM;

			out << seg3D(0,0)<< "," <<seg3D(1,0)<< "," << seg3D(2,0) << ","<<regions[eno]->regionId<<"@";
	 	}

                // need no rotation of the region (and hence all the sgments in this region) 
		else
		{
                        // translate the start point of the segment 
			seg3D << -1.0*svLength*sv[eno](0,0)+regions[eno]->midPoint(0,0),-1.0*svLength*sv[eno](1,0)+regions[eno]->midPoint(1,0),-1.0*svLength*sv[eno](2,0)+regions[eno]->zOffset;
			seg3D+=objectCM;
			out << seg3D(0,0) << "," <<seg3D(1,0)<< "," << seg3D(2,0) << ",";

			//translate the end point of the segment 
			seg3D << svLength*sv[eno](0,0)+regions[eno]->midPoint(0,0),svLength*sv[eno](1,0)+regions[eno]->midPoint(1,0),-1.0*svLength*sv[eno](2,0)+regions[eno]->zOffset;
			seg3D+=objectCM;
			out << seg3D(0,0)<< "," <<seg3D(1,0)<< "," << seg3D(2,0) << ","<<regions[eno]->regionId<<"@";
	 	}
	}

	out << endl;

	return;
}

TrajectoryVector TriMeshGeometry::extendTrajectory(Trajectory* oldtr, Direction dir)
{
	#ifdef DBG_GEOMETRY
    	cout << "DBG/GEOMETRY: Triangle::extendTrajectory() called\n";
    	#endif

        //trajectory-director
        Vector2d tDirector(0.0,0.0);

	SurfaceVector newBase(oldtr->base);

        // get copy of old trajectory region index
    	int regionIndex(oldtr->base.region->regionId),regionIndexNeigh(0);

        // forward direction: move to the end of the trajectory 
    	if(dir == ::forward)
    	{
        	newBase.x = oldtr->endPoint[1].x;
        	newBase.y = oldtr->endPoint[1].y;

                // unit vector along the direction of old trajectory
                double tDx = oldtr->endPoint[1].x-oldtr->endPoint[0].x;
    		double tDy = oldtr->endPoint[1].y-oldtr->endPoint[0].y;
    		tDirector << tDx/oldtr->length,tDy/oldtr->length;
 
                // get the neighbouring region index
        	regionIndexNeigh = oldtr->endPoint[1].nextElement;
    	}

        // backward direction: stay at the start of the trajectory
    	if(dir == backward)
    	{
                // unit vector along the direction of old trajectory
		double tDx = oldtr->endPoint[0].x-oldtr->endPoint[1].x;
    		double tDy = oldtr->endPoint[0].y-oldtr->endPoint[1].y;
                tDirector << tDx/oldtr->length,tDy/oldtr->length;

                // get the neighbouring region index
        	regionIndexNeigh = oldtr->endPoint[0].nextElement;

        	if(newBase.angle >= PI)
           	newBase.angle -= PI;
        	else
            	newBase.angle += PI;
    	}

    	// get base coordinates for a fresh trajectory (to be created)
    	callTranslator(newBase,regionIndex,regionIndexNeigh);

    	// get the new region for the fresh trajectory (to be created)
    	newBase.region = regions[regionIndexNeigh];
   	
    	double cosEdgeAngle(1.0),pCatReg(0.0);
        int s(oldtr->base.region->sideRevMap[regionIndexNeigh]);

        // cos-angle: between (encountered) triangle-edge and the trajectory-director
	double cosBendingAngle(tDirector.dot(oldtr->base.region->side[s].dir2D));

        // assign square of the cos-angle for this (encountered) triangle edge
	cosEdgeAngle = cosBendingAngle*cosBendingAngle;

        // assign pCatReg for this for this (encountered) triangle edge
	pCatReg = oldtr->base.region->side[s].pCat;
	
    	// final correction to newBase angle -> [0 ... 2*PI]
	if(newBase.angle >= 2*PI)
	newBase.angle -= 2*PI;
	else if(newBase.angle < 0)
	newBase.angle += 2*PI;

    	#ifdef DBG_GEOMETRY
    	cout << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
    	#endif

        // create a fresh trajectory
	return createAndLinkTrajectory(newBase, oldtr, dir, cosEdgeAngle, pCatReg);
}

Cartesian::Cartesian(Geometry* g, double a):Region(g, a)
{
    // necessary to use many of its member functions
    return;
}

double Cartesian::intersectionAngle(Trajectory* t1, Trajectory* t2)
{
	#ifdef DBG_ACID_TEST
	if ((t2->base.angle < 0) || (t2->base.angle > 2*PI))
		cerr << "ERROR: angle out of bounds: " << t2->base.angle << ", region=" << RegionTypeText[t2->base.region->type] << "\n";
	if ((t1->base.angle < 0) || (t1->base.angle > 2*PI))
		cerr << "ERROR: angle out of bounds: " << t1->base.angle << ", region=" << RegionTypeText[t1->base.region->type] << "\n";
	#endif

        // calculate collision angle between two interacting trajectories
	double angle(t2->base.angle - t1->base.angle);

	// angle is now in the interval -2PI..2PI
	if(angle < 0)
	angle += 2*PI; // to 0..2PI

	if(angle > PI)
	angle = 2*PI - angle;

	return (angle);
}

void Cartesian::getOrderParametersRawFlat(OrderParametersRaw& opR,vector<Vector3d>& orientation,double area)
{		
	double sin2(0.0),cos2(0.0);
	Vector3d sv(0.0,0.0,0.0);
    	double angle(0.0),length(0.0),Olength(0.0),localLength(0.);
    	double qxx(0.0),qxy(0.0),qxz(0.0),qyy(0.0),qyz(0.0),qzz(0.0);
    	double u1(0.0),u2(0.0),u3(0.0),orderLocal(0.0);
	MatrixXd QF(3,3),QC(3,3),e(3,3),O_l(2,2);

        // get the rotation matrix 
	e << 	orientation[0][0],orientation[1][0],orientation[2][0],
		orientation[0][1],orientation[1][1],orientation[2][1],
	        orientation[0][2],orientation[1][2],orientation[2][2];

    	Trajectory* tr(trajectories.first());

    	while(tr != NULL)
    	{
        	length = tr->segmentLength();

                // add up all the segmnt lengths
        	localLength += length;

        	angle = tr->base.angle;

                // scalar (local) order parameter components
        	sin2 += length*sin(2*angle);
        	cos2 += length*cos(2*angle);

        	// director for segment in local 2d-region
        	u1 = cos(angle);
        	u2 = sin(angle);
		u3 = 0.0;

        	// Q(2) tensor (local) order parameter components
        	qxx += length*(2*u1*u1);
        	qxy += length*(2*u1*u2);
        	qxz += length*(2*u1*u3);
        	qyy += length*(2*u2*u2);
        	qyz += length*(2*u2*u3);
        	qzz += length*(2*u3*u3);

        	tr = tr->next();
    	}

        // only if trace of MT segment found
	if(localLength > 0)
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
		double matrix[3][3] = {{qxx,qxy,qxz},{qxy,qyy,qyz},{qxz,qyz,qzz}};
        	double evecMat[3][3];
        	double eVal[3];
        	eigen_decomposition(matrix, evecMat, eVal);

        	int maxPos(0);

        	if(eVal[1] > eVal[0])
            	maxPos = 1;
        	if(eVal[2] > eVal[maxPos])
            	maxPos = 2;

        	for(int i=0 ; i<3 ; i++)
        	sv[i] = evecMat[i][maxPos];

		sv/= sv.norm();
		orderLocal = sqrt(sin2*sin2 + cos2*cos2)/area;

               	// global order
		qxx -= 1;
		qyy -= 1;
                QF << 	qxx,qxy,qxz,
                      	qxy,qyy,qyz,
                      	qxz,qyz,qzz;

                // make tensor transformation 
                QC = e*QF*e.transpose();
	
    		// accumulate all the Q(2) tensors components from different regions
    		opR.Qxx += localLength*QC(0,0);
    		opR.Qxy += localLength*QC(0,1);
    		opR.Qxz += localLength*QC(0,2);
    		opR.Qyy += localLength*QC(1,1);
    		opR.Qyz += localLength*QC(1,2);
    		opR.Qzz += localLength*QC(2,2);	

		// add up length of the segments from all the regions
		opR.localL += localLength;	
	}

        opR.localOrder = orderLocal;
	opR.Sv = sv;
	return;
}

void Cartesian::outputSnapshot(ostream& out)
{
	Cartesian::outputSnapshotOffset(out,0.0,0.0);
	return;
}

void Cartesian::outputOrderHeatMap(ostream& out,vector<double>& localOrder,vector<Vector3d>& sv)
{
	return;
}

void Cartesian::outputSnapshotOffset(ostream& out,double xOffset,double yOffset)
{
        Vector3d seg3D;

        // get the region 3D rotation axis 
	Vector3d cross(Q.x(),Q.y(),Q.z());
     	double crossMagnitude = cross.norm();  
  
	Quaternion<double> qr(0.0,0.0,0.0,0.0),QR(0.0,0.0,0.0,0.0);

	Trajectory* tr;
	SurfaceVector svec;
	list<Segment*>::iterator seg;

	tr = trajectories.first();

	while(tr != NULL)
	{	
		svec = tr->base;

                // record segments in the snapshots	
		seg = tr->segments.begin();
		while (seg != tr->segments.end())
		{
			svec = tr->base;

                        // move to the starting end of the segment 
			translateVector(svec, (**seg).start);

			// for backward direction retrieve actual base angle (during nucleation)
			if ((**seg).dir == backward)
			{
				svec.angle -= PI;
				if (svec.angle < 0)
				svec.angle += 2*PI;
			}

                        // region (and hence all the sgments in this region) need 3D rotation + translation
     			if(crossMagnitude>0.0)
     			{
                                // translate and rotate the start point of the segment 
        			QR.w()= 0;
				QR.x() = svec.x+midPoint(0,0);
				QR.y() = svec.y+midPoint(1,0);
				QR.z() = zOffset;

        			qr = Q.inverse() * QR * Q;
				seg3D << qr.x(),qr.y(),qr.z();
				seg3D+=geometry->objectCM;

				out << seg3D(0,0) << "," <<seg3D(1,0)<< "," << seg3D(2,0) << ",";

                                // translate and rotate the end point of the segment
				QR.w()= 0;
				QR.x() = svec.x+midPoint(0,0)+(**seg).length()*cos(svec.angle);
				QR.y() = svec.y+midPoint(1,0)+(**seg).length()*sin(svec.angle);
				QR.z() = zOffset;

        			qr = Q.inverse() * QR * Q;
        			seg3D << qr.x(),qr.y(),qr.z();
				seg3D+=geometry->objectCM;

				out << seg3D(0,0)<< "," <<seg3D(1,0)<< "," << seg3D(2,0) << ","<<regionId<<"@";
    		 	}

                        // need no rotation of the region (and hence all the sgments in this region) 
   			else
    			{
                                // translate the start point of the segment 
				seg3D << svec.x+midPoint(0,0),svec.y+midPoint(1,0),zOffset;
				seg3D+=geometry->objectCM;
				out << seg3D(0,0) << "," <<seg3D(1,0)<< "," << seg3D(2,0) << ",";
	
				//translate the end point of the segment 
        			seg3D << svec.x+midPoint(0,0)+(**seg).length()*cos(svec.angle),svec.y+midPoint(1,0)+(**seg).length()*sin(svec.angle),zOffset;
				seg3D+=geometry->objectCM;
				out << seg3D(0,0)<< "," <<seg3D(1,0)<< "," << seg3D(2,0) << ","<<regionId<<"@";
   		 	}
			seg++;
		}

		tr = tr->next();
	}
	
	return;
}

void Cartesian::translateVector(SurfaceVector& sVec, const double dist)
{
        // this can be used to move a point along a trajectory (used in making snap-shots)
	sVec.x += cos(sVec.angle)*dist;
    	sVec.y += sin(sVec.angle)*dist;
    	return;
}

void Cartesian::makeIntersectionList(Trajectory* tr1)
{
	#ifdef DBG_GEOMETRY
    	cout << "DBG/GEOMETRY: Cartesian::makeIntersectionList() called\n";
    	#endif

    	double cosDiff(0.0);
    	double denominator(0.0);
    	double temp1(0.0),temp2(0.0);
    	double bp1(0.0),bp2(0.0);
    	double cos1(0.0),cos2(0.0),sin1(0.0),sin2(0.0);

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

    	while(tr2 != NULL)
    	{
                // imagine a triangle with vertices sitting at: (a) one at the trajectory intersection point and (b) other two at the base of the two trajectories
        	if(tr1 != tr2)
        	{
                        // take cos/sin of first trajectory base angle
            		cos2 = cos(tr2->base.angle);
            		sin2 = sin(tr2->base.angle);

            		cosDiff = cos1*cos2 + sin1*sin2;
            		denominator = 1 - cosDiff*cosDiff;

                        // find the intersection point distance on the trajectory
            		if(denominator > ZERO_CUTOFF)
            		{
                		denominator = 1/denominator;
                		temp1 = (tr2->base.x  - tr1->base.x)*denominator;
                		temp2 = (tr2->base.y  - tr1->base.y)*denominator;

                		bp1 = temp1*(cos1 - cosDiff*cos2) + temp2*(sin1 - cosDiff*sin2);

                                // intersetion located within the range of trajectory length
                		if(!((bp1 < 0) || (bp1 > tr1->length)))
                		{
                    			bp2 = -temp1*(cos2 - cosDiff*cos1) - temp2*(sin2 - cosDiff*sin1);
                    			#ifdef DBG_ASSERT
                    			if ((bp2 < 0) || (bp2 > tr2->length))
                    			{
                        			cout << "ERROR: DBG/extra check: incompatible results for line intersection.\n";
                        			cout << "tr1: " << bp1 << " out of " << tr1->length << ", base=" << tr1->base.x << "," << tr1->base.y << " angle=" << tr1->base.angle << "\n";
                        			cout << "tr2: " << bp2 << " out of " << tr2->length << ", base=" << tr2->base.x << "," << tr2->base.y << " angle=" << tr2->base.angle << "\n";
                        			exit(-2);
                    			}
                    			#endif

                    			// first insert 'isThat': occupancy and other are still incorrect
                    			thatItr = tr2->intersections.insert(pair<double,Intersection>(bp2,isThat));
                    			isThis.mirror = thatItr;
                    			isThis.otherTrajectory = tr2;
                    			thisItr = tr1->intersections.insert(pair<double,Intersection>(bp1,isThis));

                    			// update the reference in the other half of the intersection
                    			thatItr->second.mirror = thisItr;

                    			// call the update routine of the other trajectory, this will also correct the occupancy number
                    			tr2->newIntersection(thatItr);
                		}
            		}
        	}

        	tr2 = tr2->next();
    	}
    	#ifdef DBG_GEOMETRY
    	cout << "DBG/GEOMETRY: " << tr1->intersections.size()-1 << " intersections\n";
    	#endif

    	return;
}

Triangle::Triangle(double elementArea, Geometry* g):Cartesian(g,elementArea)
{
        // necessary to create a full geometry, made with many triangles 
	return;
}

SurfaceVector Triangle::randomSurfaceVector()
{
	SurfaceVector v;
	double r1(0.0),r2(0.0);

        // get system pointer 
    	System* s(geometry->system);

        // random number between (0,1)
    	r1 = s->randomGen.randDblExc(1.0);
    	r2 = s->randomGen.randDblExc(1.0);

        // store sqrt(r1) to use it many times, hence avoid expensive call of sqrt()
    	double sqr1 = sqrt(r1);

    	// generate uniformly and isotropically distributed nucleation points on each triangle
    	v.x = (1 - sqr1) * Vertics[0](0,0) + (sqr1 * (1 - r2)) * Vertics[1](0,0) + (sqr1 * r2) * Vertics[2](0,0);
    	v.y = (1 - sqr1) * Vertics[0](1,0) + (sqr1 * (1 - r2)) * Vertics[1](1,0) + (sqr1 * r2) * Vertics[2](1,0);
        
        // make the unifor distribution also isotropic
    	v.angle = s->randomGen.randExc(2*PI);

    	v.region = this;

    	return v;
}

void Triangle::getTrajectoryCoordinates(SurfaceVector& sVec,double& totalLength,vector<PointATedge>& endPoint,TrajectoryVector& tVec)
{
        // copy old base coordinate
	double oldX(sVec.x);
    	double oldY(sVec.y);

    	tVec.trajectory = NULL;

    	double acos(cos(sVec.angle));
    	double asin(sin(sVec.angle));

    	// direction and base angle of the trajectory
    	if(abs(acos) < ZERO_CUTOFF)
    	{
        	if (sVec.angle < PI)
        	tVec.dir = ::forward;
        	else
        	{
        	    tVec.dir = backward;
        	    sVec.angle -= PI;
        	}
    	}	

    	else if(acos<0)
    	{
        	sVec.angle += PI;
        	tVec.dir = backward;
        	if (sVec.angle>2*PI)
        	    sVec.angle-= 2*PI;
    	}
	
    	else
        tVec.dir = ::forward;

    	// use the triangle perimeter to create a trajectory, this will to make sure to have exactly two intersections 
    	Vector2d v1(sVec.x+periMeter*acos,sVec.y+periMeter*asin);
    	Vector2d v2(sVec.x-periMeter*acos,sVec.y-periMeter*asin);

    	// intersection points of the trajectory with the triangle edges
    	int encounter(0),Neno(0);
    	for(int i=0;i<3;i++)
    	{
        	if(encounter==2)
        	    break;

        	Vector2d u1;
        	u1 = Vertics[side[i].dir[0]];

        	Vector2d u2;
        	u2 = Vertics[side[i].dir[1]];

        	Vector2d u;
        	u = u2-u1;

        	Vector2d v;
        	v = v2-v1;

        	Vector2d w;
        	w = u1-v1;

        	double sval(0.0);
        	sval = (v(1,0)*w(0,0)-v(0,0)*w(1,0))/(v(0,0)*u(1,0)-v(1,0)*u(0,0));

        	if((sval>=0.0) && (sval<=1.0))
        	{
        	    encounter+=1;
        	    Neno = sideMap[i];

                    // get triangle-edge to trajectory intersection point 
        	    endPoint.push_back(PointATedge(u1(0,0)+sval*(u2(0,0)-u1(0,0)),u1(1,0)+sval*(u2(1,0)-u1(1,0)),0.0,Neno));
        	}
    	}	

    	// base coordinate (x_base < x_other) of the trajectory 
    	if(endPoint[0].x>endPoint[1].x)
    	reverse(endPoint.begin(),endPoint.end());

    	sVec.x = endPoint[0].x;
    	sVec.y = endPoint[0].y;

    	// nucleation position (distance between old base and nucleation point) on the trajectory
    	double dx(sVec.x - oldX),dy(sVec.y - oldY);
    	tVec.pos = sqrt(dx*dx + dy*dy);

        // length of the trajectory (distance between the end points of the trajectory)
    	dx = endPoint[1].x - endPoint[0].x;
    	dy = endPoint[1].y - endPoint[0].y;
    	totalLength = sqrt(dx*dx + dy*dy);

    	// eacch trajectory should have exacty two intersection point with exactly two triangle edges
    	if(encounter!=2)
    	{
        	cout << "Problem with the number of trajectory-triangle intersections (should be == 2) ---> "<< encounter<<endl;
        	exit(-1);
    	}

    	return;
}

