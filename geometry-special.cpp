/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <fstream>
#include <iostream>
/*#include <vector>
#include <queue>
#include <string>
#include <cmath>
#include <cstdlib>
*/

#pragma warning(disable:981)
#pragma warning(disable:522)
#pragma warning(disable:383)

#include "MersenneTwister.h"
#include "corticalSim.h"


/*********************** Periodic boundary conditions **************************/

Periodic::Periodic(Coord2D inSize, System* s) : 
		Geometry(g_periodic, s, inSize.x*inSize.y), 
		size(inSize)
{
	regions.push_back(new GRectangle(inSize, this, 0));
	return;
}

void Periodic::getOrderParameters(OrderParameters& op)
{
	double orientation[3][2] = {{1, 0},{0,1},{0,0}};

	OrderParametersRaw opRaw;

	static_cast<Cartesian*>(regions[0])->getOrderParametersRawFlat(opRaw, orientation);
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;

		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	op.R = opRaw.extractR(op.Rdirector);

	return;
}


double Periodic::calculateHistogram(double hist[], bool opt)
{
	return regions[0]->calculateHistogram(hist,opt);
}

void Periodic::outputSnapshot(ostream& out)
{
	out << "canvas rectangle " << size.x << " " << size.y << "\n";
	out << "base 0 0 10 1 0 0 0 1 0\n";
	regions[0]->outputSnapshot(out);

	return;
}

TrajectoryVector Periodic::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Periodic::extendTrajectory() called\n";
#endif
	SurfaceVector newBase = oldtr->base;
	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);

	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}

	RectangleSide side = (static_cast<GRectangle*>(regions[0]))->locateSide(newBase.x, newBase.y);
	
	switch(side){
	case s_top:
		newBase.y = -0.5*size.y;
		break;
	case s_left:
		newBase.x = 0.5*size.x;
		break;
	case s_bottom:
		newBase.y = 0.5*size.y;
		break;
	case s_right:
		newBase.x = -0.5*size.x;
		break;
	};


	
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
	#endif

	return createAndLinkTrajectory(newBase, oldtr, dir, 1.0);
}	



void Periodic::calculateDiscreteAngleLengths(double *al, int *an)
{
	int i;
	double angle;
	int angleIndex;
	
	Trajectory* tr;
	list<Segment*>::iterator seg;

	for (i=0; i<system->p.discreteAngleNumber; i++)
	{
		al[i] = 0.;
		an[i]=0;
	}

    tr = regions[0]->trajectories.first();
    while (tr  != NULL)
    {
        seg = tr->segments.begin();
		angle = tr->base.angle;
		if (angle >= PI)
			angle -= PI;
		if (angle >= PI)		// include it twice, because you *can* start with 2PI + epsilon
			angle -= PI;
		angleIndex = -1;
		for (i=0; i<system->p.discreteAngleNumber; i++)  // if no angle fits the given list of nucleationAngles, angleIndex = -1 -> Some error!
		{
			if (abs(angle - system->p.nucleationAngles[i]) < ZERO_CUTOFF)
				angleIndex = i;
		}
		if (angleIndex == -1)
		{
			cout << "Error in function calculateAngleLengths - no such angle found! (angle = " << angle <<")\n";
			exit(-666);
		}
        while (seg != tr->segments.end())
        {
			if ((**seg).isLastInMT()) // only count segments at the tip of MTs -> whole MT length. ( counting: s_active, s_growing_single, s_shrinking_single, s_growing_connected, s_shrinking_connected)
				an[angleIndex]++;
			al[angleIndex] += (**seg).length();
            seg++;
        }
        tr = tr->next();
    }

    for (i=0; i<system->p.discreteAngleNumber; i++)
    {
		if (an[i])
	        al[i] /= an[i] ;
    }
    return;

}


/*********************** NxN Grid boundary conditions **************************/

Grid::Grid(Coord2D inSize, int n, System* s) : 
		Geometry(g_grid, s, inSize.x*inSize.y), 
		size(inSize),
		xNumber(n)
{
	yNumber= static_cast <int> (inSize.y/(inSize.x/n - ZERO_CUTOFF));
	Coord2D subSize(inSize.x/n, inSize.y/yNumber);
	for (int i=0 ; i<n*yNumber ; i++)
		regions.push_back(new GRectangle(subSize, this, i));
	return;
}

void Grid::getOrderParameters(OrderParameters& op)
{
	double orientation[3][2] = {{1, 0},{0,1},{0,0}};

	OrderParametersRaw opRaw;

	for (int i=0 ; i < regions.size() ; i++)	
		static_cast<Cartesian*>(regions[i])->getOrderParametersRawFlat(opRaw, orientation);
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;

		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	op.R = opRaw.extractR(op.Rdirector);

	return;
}

double Grid::calculateHistogram(double hist[], bool opt)
{
	double localLength = 0.;
	for (int i=0; i<system->p.angleHistogramBins; i++)
	{
		hist[i] = 0.;
	}

	for (int i=0; i<regions.size(); i++)
	{
		localLength += regions[i]->calculateHistogram(hist, opt,true);  // true: part of sequence. 
	}
	for (int i=0; i<system->p.angleHistogramBins; i++)
	{
		hist[i] /= localLength;
	}
	return localLength;
}

void Grid::outputSnapshot(ostream& out)
{
	int i,j;
	double dx = size.x/xNumber;
	double dy = size.y/yNumber;

	double xPos, yPos;

	yPos = 0.5*size.y - 0.5*dy;
	for (i=0 ; i<yNumber ; i++)
	{
		xPos = -0.5*size.x + 0.5*dx;
		for (j=0 ; j<xNumber ; j++)
		{
			out << "canvas rectangle " << dx << " " << dy << "\n";
			out << "base " << xPos << " " << yPos << " 10 1 0 0 0 1 0\n";
			regions[i*xNumber+j]->outputSnapshot(out);
			xPos += dx;
		}
		yPos -= dy;
	}

	return;
}

TrajectoryVector Grid::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Periodic::extendTrajectory() called\n";
#endif
	SurfaceVector newBase = oldtr->base;
	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);

	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}

	RectangleSide side = (static_cast<GRectangle*>(oldtr->base.region))->locateSide(newBase.x, newBase.y);
	
	int regionIndex;
	regionIndex = oldtr->base.region->geometryRegionIndex;
	// WARNING: really ugly - improve at a later stage!!
//	for (regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);
	int regionColumn = regionIndex%xNumber;
	int regionRow = regionIndex/xNumber;

	switch(side){
	case s_top:
		newBase.y = -0.5*size.y/yNumber;
		if (regionRow == 0)
			newBase.region = regions[xNumber*(yNumber-1) + regionColumn];
		else
			newBase.region = regions[regionIndex - xNumber];
		break;
	case s_left:
		newBase.x = 0.5*size.x/xNumber;
		if (regionColumn == 0)
			newBase.region = regions[regionIndex + xNumber - 1];
		else
			newBase.region = regions[regionIndex - 1];
		break;
	case s_bottom:
		newBase.y = 0.5*size.y/yNumber;
		if (regionRow == yNumber-1)
			newBase.region = regions[regionColumn];
		else
			newBase.region = regions[regionIndex + xNumber];
		break;
	case s_right:
		newBase.x = -0.5*size.x/xNumber;
		if (regionColumn == xNumber - 1)
			newBase.region = regions[regionIndex - xNumber + 1];
		else
			newBase.region = regions[regionIndex + 1];
		break;
	};
	
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
	#endif

	return createAndLinkTrajectory(newBase, oldtr, dir,1.0);
}	



void Grid::calculateDiscreteAngleLengths(double *al, int *an)
{
	int i;
	double angle;
	int angleIndex;
	int regionIndex;

	Trajectory* tr;
	list<Segment*>::iterator seg;

	for (i=0; i<system->p.discreteAngleNumber; i++)
	{
		al[i] = 0.;
		an[i] = 0;
	}

	for (regionIndex=0 ; regionIndex<regions.size() ; regionIndex++)
	{
		tr = regions[regionIndex]->trajectories.first();
		while (tr  != NULL)
		{
			seg = tr->segments.begin();
			angle = tr->base.angle;
			if (angle >= PI)
				angle -= PI;
			if (angle >= PI)		// include it twice, because you *can* start with 2PI + epsilon
				angle -= PI;
			angleIndex = -1;
			for (i=0; i<system->p.discreteAngleNumber; i++)  // if no angle fits the given list of nucleationAngles, angleIndex = -1 -> Some error!
			{
				if (abs(angle - system->p.nucleationAngles[i]) < ZERO_CUTOFF)
					angleIndex = i;
			}
			if (angleIndex == -1)
			{
				cout << "Error in function calculateDiscreteAngleLengths - no such angle found! (angle = " << angle <<")\n";
				exit(-666);
			}
			while (seg != tr->segments.end())
			{
				if ((**seg).isLastInMT()) // only count segments at the tip of MTs -> whole MT length. ( counting: s_active, s_growing_single, s_shrinking_single, s_growing_connected, s_shrinking_connected)
					an[angleIndex]++;
				al[angleIndex] += (**seg).length();
				seg++;
			}
			tr = tr->next();
		}
	}

	for (i=0; i<system->p.discreteAngleNumber; i++)
	{
		if (an[i])
			al[i] /= an[i] ;
	}
	return;

}


/*********************** Wormhole boundary conditions **************************/

Wormhole::Wormhole(Coord2D inSize, System* s) : 
	Geometry(g_wormhole, s, inSize.x*inSize.y), 
	size(inSize)
{
	regions.push_back(new GRectangle(inSize, this, 0));
	return;
}

void Wormhole::getOrderParameters(OrderParameters& op)
{

	double orientation[3][2] = {{1, 0},{0,1},{0,0}};

	OrderParametersRaw opRaw;

	static_cast<Cartesian*>(regions[0])->getOrderParametersRawFlat(opRaw, orientation);
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;

		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	op.R = opRaw.extractR(op.Rdirector);

	return;
}

double Wormhole::calculateHistogram(double hist[], bool opt)
{
	return	regions[0]->calculateHistogram(hist, opt);

}

void Wormhole::outputSnapshot(ostream& out)
{
	out << "canvas rectangle " << size.x << " " << size.y << "\n";
	out << "base 0 0 10 1 0 0 0 1 0\n";
	regions[0]->outputSnapshot(out);
	return;
}


TrajectoryVector Wormhole::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Wormhole::extendTrajectory() called\n";
#endif
	SurfaceVector newBase = oldtr->base;
	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);

	RectangleSide side = (static_cast<GRectangle*>(regions[0]))->locateSide(newBase.x, newBase.y);

	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}

	switch(side){
		case s_top:
			newBase.y = -0.5*size.y;
			break;
		case s_bottom:
			newBase.y = 0.5*size.y;
			break;
		case s_left:
		case s_right:
			// reflection in the y axis
			newBase.angle *= -1;
			newBase.angle += PI;
			if (newBase.angle < 0)
				newBase.angle += 2*PI;
			if (newBase.y > 0)
				newBase.y -= 0.5*size.y;
			else
				newBase.y += 0.5*size.y;
			break;
	};

#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
#endif

	return createAndLinkTrajectory(newBase, oldtr, dir, 1.0);
}	


/*********************** Cylinder definitions **************************/
/*
The Cylinder is oriented along the x axis. The three regions are body, left
cap (upside down), right cap, in that order. Both caps have the
r_accounting_dontcount and r_param_modified properties set.
*/

Cylinder::Cylinder(double len, double r , System* s) : 
	Geometry(g_cylinder, s, 2*PI*r*(len + r)), 
	radius(r),
	length(len)
{
	regions.push_back(new GRectangle(Coord2D(len,2*PI*r), this, 0));
	regions.push_back(new Disc(r, this, 1, r_accounting_dontcount, r_param_modified));
	regions.push_back(new Disc(r, this, 2, r_accounting_dontcount, r_param_modified));
	return;
}

void Cylinder::getOrderParameters(OrderParameters& op)
{
	OrderParametersRaw opRaw;
	static_cast<Cartesian*>(regions[0])->getOrderParametersRawCylinder(opRaw, radius, 0,0);

	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;

	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	double orientation[3][2] = {{0,0},{0,-1},{-1,0}};
	static_cast<Cartesian*>(regions[1])->getOrderParametersRawFlat(opRaw, orientation);
	orientation[1][1] = 1;
	static_cast<Cartesian*>(regions[2])->getOrderParametersRawFlat(opRaw, orientation);
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}

	op.R = opRaw.extractR(op.Rdirector);

	return;
}


double Cylinder::calculateHistogram(double hist[], bool opt)
{
	return	regions[0]->calculateHistogram(hist, opt);

}

void Cylinder::outputSnapshot(ostream& out)
{
	out << "canvas cylinder " << length << " " << radius << "\n";
	out << "base 0 0 0 1 0 0 0 1 0\n";
	regions[0]->outputSnapshot(out);
	out << "canvas disc " << radius << "\n";		// left cap (upside down!)
	out << "base " << -0.5*length << " 0 0 0 0 -1 0 -1 0\n";
	regions[1]->outputSnapshot(out);
	out << "canvas disc " << radius << "\n";		// right cap
	out << "base " << 0.5*length << " 0 0 0 0 -1 0 1 0\n";
	regions[2]->outputSnapshot(out);
	return;
}


TrajectoryVector Cylinder::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Cylinder::extendTrajectory() called\n";
#endif

	double tempAngle;
	double cosEdgeAngle;
	SurfaceVector newBase = oldtr->base;

	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);


	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}

	if (oldtr->base.region == regions[0])
	{
		// coming from the cylindrical part
		RectangleSide side = (static_cast<GRectangle*>(regions[0]))->locateSide(newBase.x, newBase.y);
		switch(side){
		case s_top:
			newBase.y = -PI*radius;
			cosEdgeAngle = 1.0;	// smooth transition
			break;
		case s_bottom:
			newBase.y = PI*radius;
			cosEdgeAngle = 1.0; // smooth transition
			break;
		case s_left:
			// scale y coordinate to angle (range 0..2PI) and use that to calculate
			tempAngle = newBase.y/radius + PI;
			newBase.x = cos(tempAngle)*radius;
			newBase.y = sin(tempAngle)*radius;
			newBase.angle += tempAngle;
			if (newBase.angle > 2*PI)
				newBase.angle -= 2*PI;
			newBase.region = regions[1];
			cosEdgeAngle = pow(sin(oldtr->base.angle),2);
			break;
		case s_right:
			// scale inverse y coordinate to angle (range 0..2PI) and use that to calculate
			tempAngle = PI - newBase.y/radius;
			newBase.x = cos(tempAngle)*radius;
			newBase.y = sin(tempAngle)*radius;
			newBase.angle = -PI + newBase.angle + tempAngle;
			if (newBase.angle < 0)
				newBase.angle += 2*PI;
			else if (newBase.angle > 2*PI)
				newBase.angle -= 2*PI;
			newBase.region = regions[2];
			cosEdgeAngle = pow(sin(oldtr->base.angle),2);
			break;
		};

	}
	else if (oldtr->base.region == regions[1])
	{
		// left cap
		tempAngle = atan2(newBase.y, newBase.x);
		newBase.x = -0.5*length;
		if (tempAngle > 0)
			newBase.y = radius*(-PI + tempAngle);
		else
			newBase.y = radius*(PI + tempAngle);
		newBase.angle -= tempAngle;
		if (newBase.angle < 0)
			newBase.angle += 2*PI;
		else if (newBase.angle > 2*PI)
			newBase.angle -= 2*PI;
		newBase.region = regions[0];
		cosEdgeAngle = pow(sin(newBase.angle),2);
	}
	else
	{
		// right cap
		tempAngle = atan2(newBase.y, newBase.x);
		newBase.x = 0.5*length;
		if (tempAngle >= 0)
			newBase.y = radius*(PI - tempAngle);
		else
			newBase.y = -(PI + tempAngle)*radius;
		newBase.angle = PI - tempAngle + newBase.angle;
		if (newBase.angle < 0)
			newBase.angle += 2*PI;
		else if (newBase.angle > 2*PI)
			newBase.angle -= 2*PI;
		newBase.region = regions[0];
		cosEdgeAngle = pow(sin(newBase.angle),2);
	}

	// this part is generic

	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
	#endif

	return createAndLinkTrajectory(newBase, oldtr, dir, cosEdgeAngle);
}	

/*********************** GridCylinder definitions **************************/
/*
The GridCylinder is oriented along the x axis. The regions are left
cap (upside down), right cap, body rectangles, in that order. Both caps have the
r_accounting_dontcount and r_param_modified properties set.
*/

GridCylinder::GridCylinder(double len, double r, int n, System* s) : 
	Geometry(g_gridcylinder, s, 2*PI*r*(len + r)), 
	radius(r),
	length(len),
	number(n)
{
	gridLenSize = len/n;
	gridRNumber = static_cast <int> (2*PI*r/ gridLenSize);
	gridRSize = 2*PI*r / gridRNumber;

	regions.push_back(new Disc(r, this, 0, r_accounting_dontcount, r_param_modified));  //left cap
	regions.push_back(new Disc(r, this, 1, r_accounting_dontcount, r_param_modified));	//right cap
	for (int i=0; i< n*gridRNumber; i++)
		regions.push_back(new GRectangle(Coord2D(gridLenSize,gridRSize), this, i+2));
	cout << "gridcylinder length: " << len << " = " << n << " x " << gridLenSize << " = " << n * gridLenSize << "\n";
	cout << "gridcylinder circumference: " << 2*PI*r << " = " << gridRNumber << " x " << gridRSize << " = " << gridRNumber * gridRSize << "\n";
	cout << "Number of regions: " << regions.size() << "\n";
	return;
}

void GridCylinder::getOrderParameters(OrderParameters& op)
{
	OrderParametersRaw opRaw;
	for(int i=2; i < regions.size(); i++)
	{
		int regionColumn = (i-2)%number;
		int regionRow = (i-2)/number;
		static_cast<Cartesian*>(regions[i])->getOrderParametersRawCylinder(opRaw, radius, (regionColumn-(number-1)/2.)*gridLenSize, -(regionRow-(gridRNumber-1)/2.)*gridRSize);
	}
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;

	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	double orientation[3][2] = {{0,0},{0,-1},{-1,0}};
	static_cast<Cartesian*>(regions[0])->getOrderParametersRawFlat(opRaw, orientation);
	orientation[1][1] = 1;
	static_cast<Cartesian*>(regions[1])->getOrderParametersRawFlat(opRaw, orientation);
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}

	op.R = opRaw.extractR(op.Rdirector);

	return;
}

double GridCylinder::calculateHistogram(double hist[], bool opt)
{
	double localLength = 0.;
	for (int i=0; i<system->p.angleHistogramBins; i++)
	{
		hist[i] = 0.;
	}
	for (int i=2; i<regions.size(); i++)
	{
		localLength += regions[i]->calculateHistogram(hist, opt, true);  // true: part of sequence. 
	}
	for (int i=0; i<system->p.angleHistogramBins; i++)
	{
		hist[i] /= localLength;
	}
	return localLength;
}

void GridCylinder::outputSnapshot(ostream& out)
{
	out << "canvas disc " << radius << "\n";		// left cap (upside down!)
	out << "base " << -0.5*length << " 0 0 0 0 -1 0 -1 0\n";
	regions[0]->outputSnapshot(out);
	out << "canvas disc " << radius << "\n";		// right cap
	out << "base " << 0.5*length << " 0 0 0 0 -1 0 1 0\n";
	regions[1]->outputSnapshot(out);
	out << "canvas cylinder " << length << " " << radius << "\n";
	out << "base 0 0 0 1 0 0 0 1 0\n";
	for(int i=2; i < regions.size(); i++)
	{
		int regionColumn = (i-2)%number;
		int regionRow = (i-2)/number;
		static_cast<Cartesian*>(regions[i])->outputSnapshotOffset(out, (regionColumn-(number-1)/2.)*gridLenSize, -(regionRow-(gridRNumber-1)/2.)*gridRSize);
	}
	return;
}


TrajectoryVector GridCylinder::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: GridCylinder::extendTrajectory() called\n";
#endif

	double tempAngle;
	double cosEdgeAngle;
	SurfaceVector newBase = oldtr->base;
	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);


	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}

	if (oldtr->base.region == regions[0])
	{
		// left cap
		double tempY;
		int regionIndex;
		tempAngle = atan2(newBase.y, newBase.x);
		newBase.x = -0.5*gridLenSize;
		if (tempAngle > 0)
			tempY = radius*(-PI + tempAngle);
		else
			tempY = radius*(PI + tempAngle);
		newBase.angle -= tempAngle;
		if (newBase.angle < 0)
			newBase.angle += 2*PI;
		else if (newBase.angle > 2*PI)
			newBase.angle -= 2*PI;
		newBase.y = fmod(tempY + 3*PI*radius,gridRSize) - gridRSize * 0.5;
		regionIndex = 2 + number * static_cast <int> ((PI*radius - tempY )/ gridRSize); // columnNumber: 0.
		newBase.region = regions[regionIndex];
		cosEdgeAngle = pow(sin(newBase.angle),2);
	}
	else if (oldtr->base.region == regions[1]) 
	{
		// right cap
		double tempY;
		int regionIndex;
		tempAngle = atan2(newBase.y, newBase.x);
		newBase.x = 0.5*gridLenSize;
		if (tempAngle >= 0)
			tempY = radius*(PI - tempAngle);
		else
			tempY = -(PI + tempAngle)*radius;
		newBase.angle = PI - tempAngle + newBase.angle;
		if (newBase.angle < 0)
			newBase.angle += 2*PI;
		else if (newBase.angle > 2*PI)
			newBase.angle -= 2*PI;
		newBase.y = fmod(tempY + 3*PI*radius,gridRSize) - gridRSize * 0.5;
		regionIndex = 1 + number + number * static_cast <int> ((PI*radius - tempY )/ gridRSize); // columnNumber: number - 1. 
		newBase.region = regions[regionIndex];
		cosEdgeAngle = pow(sin(newBase.angle),2);
	}
	else
	{
		// coming from the cylindrical part
		int regionIndex;
		// WARNING: really ugly - improve at a later stage!!
		// requires a reverse lookup: *Region -> regionIndex
		// perhaps store the regionindex inside the region?
//		for (regionIndex=2 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);
		regionIndex = oldtr->base.region->geometryRegionIndex;

		int regionColumn = (regionIndex-2)%number;
		int regionRow = (regionIndex-2)/number;
		RectangleSide side = (static_cast<GRectangle*>(regions[regionIndex]))->locateSide(newBase.x, newBase.y);
		switch(side){
			case s_top:
				newBase.y = -0.5*gridRSize;
				cosEdgeAngle = 1.0;
				if (regionRow == 0)
					newBase.region = regions[number*(gridRNumber-1) + regionColumn + 2];
				else
					newBase.region = regions[regionIndex - number];
				break;
			case s_bottom:
				newBase.y = 0.5*gridRSize;
				cosEdgeAngle = 1.0;
				if (regionRow == gridRNumber-1)
					newBase.region = regions[regionColumn + 2 ];
				else
					newBase.region = regions[regionIndex + number];
				break;
			case s_left:
				if (regionColumn > 0)
				{
					newBase.x = 0.5*gridLenSize;
					newBase.region = regions[regionIndex - 1];
					cosEdgeAngle = 1.0;
				}
				else
				{
					// scale y coordinate to angle (range 0..2PI) and use that to calculate
					tempAngle = (newBase.y -(regionRow-(gridRNumber-1)/2.)*gridRSize)/radius + PI ; 
					newBase.x = cos(tempAngle)*radius;
					newBase.y = sin(tempAngle)*radius;
					newBase.angle += tempAngle;
					if (newBase.angle > 2*PI)
						newBase.angle -= 2*PI;
					newBase.region = regions[0];
					cosEdgeAngle = pow(sin(oldtr->base.angle),2);
				}
				break;
			case s_right:
				if (regionColumn < number - 1)
				{
					newBase.x = -0.5*gridLenSize;
					newBase.region = regions[regionIndex + 1];
					cosEdgeAngle = 1.0;
				}
				else
				{
					// scale inverse y coordinate to angle (range 0..2PI) and use that to calculate
					tempAngle = PI - (newBase.y -(regionRow-(gridRNumber-1)/2.)*gridRSize)/radius; 
					newBase.x = cos(tempAngle)*radius;
					newBase.y = sin(tempAngle)*radius;
					newBase.angle = -PI + newBase.angle + tempAngle;
					if (newBase.angle < 0)
						newBase.angle += 2*PI;
					else if (newBase.angle > 2*PI)
						newBase.angle -= 2*PI;
					newBase.region = regions[1];
					cosEdgeAngle = pow(sin(oldtr->base.angle),2);
				}
				break;
		};

	}

	// this part is generic

#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
#endif

	return createAndLinkTrajectory(newBase, oldtr, dir, cosEdgeAngle);
}	


/*********************** Pancake boundary conditions **************************/

Pancake::Pancake(double r , System* s) : 
		Geometry(g_pancake, s, 2*PI*r*r), 
		radius(r)
{
	regions.push_back(new Disc(r, this, 0));
	return;
}

void Pancake::getOrderParameters(OrderParameters& op)
{
	double orientation[3][2] = {{1, 0},{0,1},{0,0}};

	OrderParametersRaw opRaw;

	static_cast<Cartesian*>(regions[0])->getOrderParametersRawFlat(opRaw, orientation);
	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;

		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	op.R = opRaw.extractR(op.Rdirector);

	return;
}


double Pancake::calculateHistogram(double hist[], bool opt)
{
	return	regions[0]->calculateHistogram(hist, opt);

}

void Pancake::outputSnapshot(ostream& out)
{
	out << "canvas disc " << radius << "\n";
	out << "base 0 0 0 1 0 0 0 1 0\n";
	regions[0]->outputSnapshot(out);
	return;
}


TrajectoryVector Pancake::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Cylinder::extendTrajectory() called\n";
#endif

	SurfaceVector newBase = oldtr->base;
	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);

	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}
	
	newBase.x *= -1;
	newBase.y *= -1;


	// this part is generic

	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
	#endif

	return createAndLinkTrajectory(newBase, oldtr, dir, 1.0);
}	


/*********************** Box **************************/
/*
 * 
 * 
 */ 
 
Box::Box(double xin, double yin, double zin, System* s) 
	: Geometry(g_box, s, 2*(xin*yin + xin*zin + yin*zin)), 
	  x(xin), 
	  y(yin), 
	  z(zin)
{
	regions.push_back(new GRectangle(Coord2D(x,z), this, 0, r_accounting_dontcount, r_param_modified));	// top
	regions.push_back(new GRectangle(Coord2D(x,y), this, 1));	// front
	regions.push_back(new GRectangle(Coord2D(z,y), this, 2));	// right
	regions.push_back(new GRectangle(Coord2D(x,y), this, 3));	// back
	regions.push_back(new GRectangle(Coord2D(z,y), this, 4));	// left
	regions.push_back(new GRectangle(Coord2D(x,z), this, 5, r_accounting_dontcount, r_param_modified));	// bottom
	return;
}



void Box::getOrderParameters(OrderParameters& op)
{

	OrderParametersRaw opRaw;

	// first, collect data on the surrounding faces (excluding top and bottom)
	// order: front, right, back, left
	double orientation[3][2] = {{1, 0},{0,1},{0,0}};
	static_cast<Cartesian*>(regions[1])->getOrderParametersRawFlat(opRaw, orientation);
	orientation[0][0] = 0;
	orientation[2][0] = -1;
	static_cast<Cartesian*>(regions[2])->getOrderParametersRawFlat(opRaw, orientation);
	orientation[2][0] = 0;
	orientation[0][0] = -1;
	static_cast<Cartesian*>(regions[3])->getOrderParametersRawFlat(opRaw, orientation);
	orientation[0][0] = 0;
	orientation[2][0] = 1;
	static_cast<Cartesian*>(regions[4])->getOrderParametersRawFlat(opRaw, orientation);

	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.si2 /= opRaw.localL;
		opRaw.si4 /= opRaw.localL;
		opRaw.co2 /= opRaw.localL;
		opRaw.co4 /= opRaw.localL;
		opRaw.si2Opt /= opRaw.localLOpt;
		opRaw.si4Opt /= opRaw.localLOpt;
		opRaw.co2Opt /= opRaw.localLOpt;
		opRaw.co4Opt /= opRaw.localLOpt;
	}
	op.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
	op.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
	op.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
	op.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
	op.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
	op.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
	op.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
	op.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );

	orientation[1][1] = 0;
	orientation[2][1] = -1;
	orientation[2][0] = 0;
	orientation[0][0] = 1;
	static_cast<Cartesian*>(regions[0])->getOrderParametersRawFlat(opRaw, orientation);
	orientation[2][1] = 1;
	static_cast<Cartesian*>(regions[5])->getOrderParametersRawFlat(opRaw, orientation);

	if (opRaw.localL > ZERO_CUTOFF)
	{
		opRaw.Qxx /= opRaw.localL;
		opRaw.Qxy /= opRaw.localL;
		opRaw.Qxz /= opRaw.localL;
		opRaw.Qyy /= opRaw.localL;
		opRaw.Qyz /= opRaw.localL;
		opRaw.Qzz /= opRaw.localL;

		opRaw.isoWeights[0] /= area;
		opRaw.isoWeights[1] /= area;
		opRaw.isoWeights[2] /= area;
	}


	op.R = opRaw.extractR(op.Rdirector);

	return;
}


double Box::calculateHistogram(double hist[], bool opt)
{
	return	regions[1]->calculateHistogram(hist, opt);
	
}



void Box::outputSnapshot(ostream& out)
{
	out << "canvas rectangle " << x << " " << z << "\n";	// top
	out << "base 0 " << 0.5*y << " 0 1 0 0 0 0 -1\n";
	regions[0]->outputSnapshot(out);
	out << "canvas rectangle " << x << " " << y << "\n";	// front
	out << "base 0 0 " << 0.5*z << " 1 0 0 0 1 0\n";
	regions[1]->outputSnapshot(out);
	out << "canvas rectangle " << z << " " << y << "\n";	// right
	out << "base " << 0.5*x << " 0 0 0 0 -1 0 1 0\n";
	regions[2]->outputSnapshot(out);
	out << "canvas rectangle " << x << " " << y << "\n";	// back
	out << "base 0 0 " << -0.5*z << " -1 0 0 0 1 0\n";
	regions[3]->outputSnapshot(out);
	out << "canvas rectangle " << z << " " << y << "\n";	// left
	out << "base " << -0.5*x << " 0 0 0 0 1 0 1 0\n";
	regions[4]->outputSnapshot(out);
	out << "canvas rectangle " << x << " " << z << "\n";	// bottom
	out << "base 0 " << -0.5*y << " 0 1 0 0 0 0 1\n";
	regions[5]->outputSnapshot(out);
	return;
}


TrajectoryVector Box::extendTrajectory(Trajectory* oldtr, Direction dir)
{
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Cylinder::extendTrajectory() called\n";
	#endif

	double cosEdgeAngle;

	SurfaceVector newBase = oldtr->base;
	if (dir == ::forward)
		oldtr->base.region->translateVector(newBase,oldtr->length);

	if (dir == backward)
	{
		if (newBase.angle >= PI)
			newBase.angle -= PI;
		else
			newBase.angle += PI;
	}
	
	RectangleSide side = (static_cast<GRectangle*>(oldtr->base.region))->locateSide(newBase.x, newBase.y);
	if (oldtr->base.region == regions[0])
	{
		// coming from the top
		switch(side){
		case s_top:
			// go to back
			newBase.x *= -1;
			newBase.y = 0.5*y;
			newBase.angle += PI;
			newBase.region = regions[3];
			break;
		case s_bottom:
			// go to front
			newBase.y = 0.5*y;
			newBase.region = regions[1];
			break;
		case s_left:
			// go to left
			newBase.x = -newBase.y;
			newBase.y = 0.5*y;
			newBase.angle += 0.5*PI;
			newBase.region = regions[4];
			break;
		case s_right:
			// go to right
			newBase.x = newBase.y;
			newBase.y = 0.5*y;
			newBase.angle -= 0.5*PI;
			newBase.region = regions[2];
			break;
		};
	}
	else if (oldtr->base.region == regions[1])
	{
		// front
		switch(side){
		case s_top:
			// go to top
			newBase.y = -0.5*z;
			newBase.region = regions[0];
			break;
		case s_bottom:
			// go to bottom
			newBase.y = 0.5*z;
			newBase.region = regions[5];
			break;
		case s_left:
			// go to left
			newBase.x = 0.5*z;
			newBase.region = regions[4];
			break;
		case s_right:
			// go to right
			newBase.x = -0.5*z;
			newBase.region = regions[2];
			break;
		};
	}
	else if (oldtr->base.region == regions[2])
	{
		// coming from the right...
		switch(side){
		case s_top:
			// to the top
			newBase.y = newBase.x;
			newBase.x = 0.5*x;
			newBase.angle += 0.5*PI;
			newBase.region = regions[0];
			break;
		case s_bottom:
			// to the bottom
			newBase.y = -newBase.x;
			newBase.x = 0.5*x;
			newBase.angle -= 0.5*PI;
			newBase.region = regions[5];
			break;
		case s_left:
			// to the front
			newBase.x = 0.5*x;
			newBase.region = regions[1];
			break;
		case s_right:
			// to the back
			newBase.x = -0.5*x;
			newBase.region = regions[3];
			break;
		};
	}
	else if (oldtr->base.region == regions[3])
	{
		// coming from the back...
		switch(side){
		case s_top:
			// to the top
			newBase.x *= -1;
			newBase.y = 0.5*z;
			newBase.angle += PI;
			newBase.region = regions[0];
			break;
		case s_bottom:
			// to the bottom
			newBase.x *= -1;
			newBase.y = -0.5*z;
			newBase.angle += PI;
			newBase.region = regions[5];
			break;
		case s_left:
			// to the right
			newBase.x = 0.5*z;
			newBase.region = regions[2];
			break;
		case s_right:
			// to the left
			newBase.x = -0.5*z;
			newBase.region = regions[4];
			break;
		};
	}
	else if (oldtr->base.region == regions[4])
	{
		// coming from the left...
		switch(side){
		case s_top:
			// to the top
			newBase.y = -newBase.x;
			newBase.x = -0.5*x;
			newBase.angle -= 0.5*PI;
			newBase.region = regions[0];
			break;
		case s_bottom:
			// to the bottom
			newBase.y = newBase.x;
			newBase.x = -0.5*x;
			newBase.angle += 0.5*PI;
			newBase.region = regions[5];
			break;
		case s_left:
			// to the back
			newBase.x = 0.5*x;
			newBase.region = regions[3];
			break;
		case s_right:
			// to the front
			newBase.x = -0.5*x;
			newBase.region = regions[1];
			break;
		};
	}
	else // if (oldtr->base.region == regions[5])
	{
		// coming from the bottom...
		switch(side){
		case s_top:
			// to the front
			newBase.y = -0.5*y;
			newBase.region = regions[1];
			break;
		case s_bottom:
			// to the back
			newBase.x *= -1;
			newBase.y = -0.5*y;
			newBase.angle += PI;
			newBase.region = regions[3];
			break;
		case s_left:
			// to the left
			newBase.x = newBase.y;
			newBase.y = -0.5*y;
			newBase.angle -= 0.5*PI;
			newBase.region = regions[4];
			break;
		case s_right:
			// to the right
			newBase.x = -newBase.y;
			newBase.y = -0.5*y;
			newBase.angle += 0.5*PI;
			newBase.region = regions[2];
			break;
		};
	}

	switch(side){
	case s_top:
	case s_bottom:
		cosEdgeAngle = pow(cos(oldtr->base.angle),2);
		break;
	case s_left:
	case s_right:
		cosEdgeAngle = pow(sin(oldtr->base.angle),2);
		break;
	};
	
	if (newBase.angle >= 2*PI)
		newBase.angle -= 2*PI;
	else if (newBase.angle < 0)
		newBase.angle += 2*PI;

	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Extended trajectory from region type " << RegionTypeText[oldtr->base.region->type] << " to " << RegionTypeText[newBase.region->type] << "\n";
	#endif

	return createAndLinkTrajectory(newBase, oldtr, dir, cosEdgeAngle);
}	



/********************************** CARTESIAN REGION FUNCTIONS ****************************/
Cartesian::Cartesian(RegionType r, Geometry* g, int idx, double a, RegionAccountingType acc, RegionParameterType par) : 
		Region(r, g, idx, a, acc, par)
{
	return;
}

double Cartesian::intersectionAngle(Trajectory* t1, Trajectory* t2)
/* output is expected to lie in the interval 0..PI
 * 
 */
{
	#ifdef DBG_ACID_TEST
	if ((t2->base.angle < 0) || (t2->base.angle > 2*PI))
		cerr << "ERROR: angle out of bounds: " << t2->base.angle << ", region=" << RegionTypeText[t2->base.region->type] << "\n";
	if ((t1->base.angle < 0) || (t1->base.angle > 2*PI))
		cerr << "ERROR: angle out of bounds: " << t1->base.angle << ", region=" << RegionTypeText[t1->base.region->type] << "\n";
	
	#endif
	
	double angle = t2->base.angle - t1->base.angle;
	// angle is now in the interval -2PI..2PI
	if (angle < 0)
		angle += 2*PI; // to 0..2PI
	
	if (angle > PI)
		angle = 2*PI - angle;
	
	return angle;	
}

void Cartesian::getOrderParametersRawFlat(OrderParametersRaw& opR, double orientation[3][2])
{
		int i;

        double sin2,cos2,cos4,sin4,Osin2,Ocos2,Ocos4,Osin4;
        double angle;
        double length, Olength;
        double localLength, OlocalLength;

		double u1,u2;
		double v[3];
		double qxx,qxy,qxz,qyy,qyz,qzz;
		qxx=qxy=qxz=qyy=qyz=qzz=0;


        Trajectory* tr;

        sin2=cos2=sin4=cos4=0;
        Osin2=Ocos2=Osin4=Ocos4=0;
        localLength = 0.;
        OlocalLength = 0.;
        tr = trajectories.first();



        while (tr != NULL)
        {
                length = tr->segmentLength();
                localLength += length;
                Olength = tr->coveredLength();
                OlocalLength += Olength;
                angle = tr->base.angle;
                sin2 += length*sin(2*angle);
                cos2 += length*cos(2*angle);
                sin4 += length*sin(4*angle);
                cos4 += length*cos(4*angle);
                Osin2 += Olength*sin(2*angle);
                Ocos2 += Olength*cos(2*angle);
                Osin4 += Olength*sin(4*angle);
                Ocos4 += Olength*cos(4*angle);

				u1 = cos(angle);
				u2 = sin(angle);
				for (i=0 ; i<3 ; i++)
					v[i] = orientation[i][0]*u1 + orientation[i][1]*u2;

				qxx += length*v[0]*v[0];
				qxy += length*v[0]*v[1];
				qxz += length*v[0]*v[2];
				qyy += length*v[1]*v[1];
				qyz += length*v[1]*v[2];
				qzz += length*v[2]*v[2];

				tr = tr->next();

        }
        opR.si2 += sin2;
        opR.si4 += sin4;
        opR.co2 += cos2;
        opR.co4 += cos4;
		opR.localL += localLength;
        opR.si2Opt += Osin2;
        opR.si4Opt += Osin4;
        opR.co2Opt += Ocos2;
        opR.co4Opt += Ocos4;
		opR.localLOpt += OlocalLength;

		opR.Qxx += qxx;
		opR.Qxy += qxy;
		opR.Qxz += qxz;
		opR.Qyy += qyy;
		opR.Qyz += qyz;
		opR.Qzz += qzz;

		for (i=0 ; i<3 ; i++)
		{
			opR.isoWeights[i] += area*0.5*sqrt(pow(orientation[i][0],2) + pow(orientation[i][1],2));
		}

        return;
}

void Cartesian::getOrderParametersRawCylinder(OrderParametersRaw& opR, double radius, double xOffset, double yOffset)
/* WARNING: the calculations of Qab and the isoWeights assume the default orientation of the cylinder!
*/
{

        double sin2,cos2,cos4,sin4,Osin2,Ocos2,Ocos4,Osin4;
        double angle;
        double length, Olength;
        double localLength, OlocalLength;

		double baseAngle;
		double qxx,qxy,qxz,qyy,qyz,qzz;
		qxx=qxy=qxz=qyy=qyz=qzz=0;
		list<Segment*>::iterator seg;


        Trajectory* tr;

        sin2=cos2=sin4=cos4=0;
        Osin2=Ocos2=Osin4=Ocos4=0;
        localLength = 0.;
        OlocalLength = 0.;
        tr = trajectories.first();



        while (tr != NULL)
        {
                length = tr->segmentLength();
                localLength += length;
                Olength = tr->coveredLength();
                OlocalLength += Olength;
                angle = tr->base.angle;
                sin2 += length*sin(2*angle);
                cos2 += length*cos(2*angle);
                sin4 += length*sin(4*angle);
                cos4 += length*cos(4*angle);
                Osin2 += Olength*sin(2*angle);
                Ocos2 += Olength*cos(2*angle);
                Osin4 += Olength*sin(4*angle);
                Ocos4 += Olength*cos(4*angle);

	
				seg = tr->segments.begin();
				while (seg != tr->segments.end())
				{
					length = (**seg).length();
					if ((*seg)->dir == ::forward)
						baseAngle = (yOffset + tr->base.y + sin(angle)*((*seg)->start))/radius;
					else
						baseAngle = (yOffset + tr->base.y + sin(angle)*((*seg)->end))/radius;

					qxx += length*pow(cos(angle),2);
					qxy += radius*cos(angle)*(- sin(baseAngle) + sin(baseAngle+length*sin(angle)/radius));
					qxz += radius*cos(angle)*(- cos(baseAngle) + cos(baseAngle+length*sin(angle)/radius));
					qyy += 0.25*sin(angle)*(2*length*sin(angle) + radius*(-sin(2*baseAngle)+sin(2*(baseAngle+length*sin(angle)/radius))));
					qyz += 0.5*radius*sin(angle)*(-pow(cos(baseAngle),2) +pow(cos(baseAngle+length*sin(angle)/radius),2));
					qzz += 0.25*sin(angle)*(2*length*sin(angle) + radius*(sin(2*baseAngle)-sin(2*(baseAngle+length*sin(angle)/radius)))); 

					seg++;
				}

				tr = tr->next();

        }
        opR.si2 += sin2;
        opR.si4 += sin4;
        opR.co2 += cos2;
        opR.co4 += cos4;
		opR.localL += localLength;
        opR.si2Opt += Osin2;
        opR.si4Opt += Osin4;
        opR.co2Opt += Ocos2;
        opR.co4Opt += Ocos4;
		opR.localLOpt += OlocalLength;

		opR.Qxx += qxx;
		opR.Qxy += qxy;
		opR.Qxz += qxz;
		opR.Qyy += qyy;
		opR.Qyz += qyz;
		opR.Qzz += qzz;

		opR.isoWeights[0] += area*0.5;
		opR.isoWeights[1] += area*0.25;
		opR.isoWeights[2] += area*0.25;

        return;
}


double Cartesian::calculateHistogram(double hist[], bool optical, bool partOfSequence)
{
	int i;
	double angle;
	double localLength, lengthInc;

	Trajectory* tr;
	int bin;
	const double binsize=PI/geometry->system->p.angleHistogramBins;		// vanaf 7- 4: hoeken over -1/2 PI - 1/2 PI ( 0 in het midden! )
	// const double binsize = 2.*PI/angleHistogram;

	localLength = 0.;
	if (! partOfSequence)
	{
		for (i=0; i<geometry->system->p.angleHistogramBins; i++)
		{
			hist[i] = 0.;
		}
	}

	tr = trajectories.first();
	while (tr != NULL)
	{

		angle = tr->base.angle;
		bin = (static_cast <int> (floor((angle+2.0*PI)/binsize + 0.5))%geometry->system->p.angleHistogramBins);
#ifdef DBG_ASSERT
		if (bin < 0 || bin >= geometry->system->p.angleHistogramBins)
		{
			cerr << "ERROR: bin out of range. [angle=" << angle << ", bin=" << bin << "]. Exiting.\n";
			exit(-1);
		}
#endif
    if (optical)
    {
      lengthInc = tr->coveredLength();
    }
    else
    {
      lengthInc = tr->segmentLength();
    }

		hist[bin] += lengthInc;
		localLength += lengthInc;
		tr = tr->next();
	}
	if (! partOfSequence)
	{
		if (localLength > ZERO_CUTOFF)		
		{
			for (i=0; i<geometry->system->p.angleHistogramBins; i++)
			{
				hist[i] /= localLength;
			}
		}
	}
	return localLength;
}

void Cartesian::outputSnapshot(ostream& out)
{
	Cartesian::outputSnapshotOffset(out, 0., 0.);
	return;
}
void Cartesian::outputSnapshotOffset(ostream& out, double xOffset, double yOffset)
{
	Trajectory* tr;
	list<Segment*>::iterator seg;
	SurfaceVector svec;

	tr = trajectories.first();
	while (tr != NULL)
	{
		svec = tr->base;
		// uncomment the following code to record trajectories in the snapshots
		//		out << "l " << graphics_trajectory << " ";
		//		out << svec.x << " " << svec.y << " " << svec.angle << " " << tr->length;
		//		out << "\n";

		seg = tr->segments.begin();
		while (seg != tr->segments.end())
		{
			svec = tr->base;
			svec.x += xOffset;
			svec.y += yOffset;
			translateVector(svec, (**seg).start);
			if ((**seg).dir == backward)
			{
				svec.angle -= PI;
				if (svec.angle < 0)
					svec.angle += 2*PI;
			}
			out << "l " << graphics_mt << " ";
			out << svec.x << " " << svec.y << " " << svec.angle << " " << (**seg).length();
			out << "\n";

			if ((**seg).isFirstInMT())
			{
				out << "p " << graphics_minus << " ";
				out << svec.x << " " << svec.y ;
				out << "\n";
			}
			if ((**seg).isLastInMT())
			{
				translateVector(svec, (**seg).length());
				out << "p " << graphics_plus << " ";
				out << svec.x << " " << svec.y ;
				out << "\n";
			}

			seg++;
		}
		tr = tr->next();
	}

	return;
}



void Cartesian::translateVector(SurfaceVector& sVec, const double dist)
{
	sVec.x += cos(sVec.angle)*dist;
	sVec.y += sin(sVec.angle)*dist;
	
	return;	
}

void Cartesian::makeIntersectionList(Trajectory* tr1)
{
#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Cartesian::makeIntersectionList() called\n";
#endif

	Trajectory* tr2;
	
	double cosDiff;
	double denominator;
	double temp1, temp2;
	double bp1,bp2;
	double cos1,cos2,sin1,sin2;

	cos1 = cos(tr1->base.angle);
	sin1 = sin(tr1->base.angle);

	Intersection isThis;
	Intersection isThat;
	IntersectionItr thisItr;
	IntersectionItr thatItr;
	
	isThis.occupancy = 0;
	isThat.occupancy = 0;
	isThat.otherTrajectory = tr1;

	tr2 = trajectories.first();
	while (tr2 != NULL)
	{
		if (tr1 != tr2)
		{
			cos2 = cos(tr2->base.angle);
			sin2 = sin(tr2->base.angle);
			
			cosDiff = cos1*cos2 + sin1*sin2;
			denominator = 1 - cosDiff*cosDiff;
		
			if (denominator > ZERO_CUTOFF)
			{
				denominator = 1/denominator;
				temp1 = (tr2->base.x  - tr1->base.x)*denominator;
				temp2 = (tr2->base.y  - tr1->base.y)*denominator;
				bp1 = temp1*(cos1 - cosDiff*cos2) + temp2*(sin1 - cosDiff*sin2);
							
				if (!((bp1 < 0) || (bp1 > tr1->length)))
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
								
					// now update the reference in the other half of the intersection
					thatItr->second.mirror = thisItr;
					// call the update routine of the other trajectory. This will also correct the occupancy number
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



/*********************************** RECTANGULAR REGION FUNCTIONS *******************************/

GRectangle::GRectangle(Coord2D inSize, Geometry* g, int idx, RegionAccountingType acc, RegionParameterType par) : 
		Cartesian(r_rectangle, g, idx, inSize.x*inSize.y, acc, par), 
		size(inSize)
{
	return;
}


SurfaceVector GRectangle::randomSurfaceVector()
{
	SurfaceVector v;
	System* s = geometry->system;
	// make sure that no microtubules are nucleated exactly on a region boundary
	// this can lead to event ordering problems
	v.x = s->randomGen.randDblExc(size.x) - 0.5*size.x;
	v.y = s->randomGen.randDblExc(size.y) - 0.5*size.y;
	v.angle = s->randomGen.randExc(2*PI);
	v.region = this;

	return v;	
}


void GRectangle::getTrajectoryCoordinates(SurfaceVector& sVec, double& totalLength, TrajectoryVector& tVec)
{
	#ifdef DBG_ACID_TEST
	if ((abs(sVec.x) > (0.5 + ZERO_CUTOFF)*size.x) || (abs(sVec.y) > (0.5 + ZERO_CUTOFF)*size.y))
	{
		cerr << "GRectangle::getTrajectoryCoordinates() received an out of bounds coordinate.\n";
		exit(-1);
	}
	#endif

	tVec.trajectory = NULL;
	double acos = cos(sVec.angle);
	double asin = sin(sVec.angle);


	// prevent divide by zeroes
	if (abs(acos) < ZERO_CUTOFF)
	{
		// this one is a little tricky, because it must be done consistently with the direction
		tVec.pos = sVec.y + 0.5*size.y;
		sVec.y = -0.5*size.y;
		totalLength = size.y;
		if (sVec.angle < PI)
			tVec.dir = ::forward;
		else
		{
			tVec.dir = backward;
			sVec.angle -= PI;
		}
		return;
	}

	if (acos<0)
	{
		acos *= -1;
		asin *=-1;
		sVec.angle += PI;
		tVec.dir = backward;
		if (sVec.angle>2*PI)
			sVec.angle-= 2*PI;
	}
	else
		tVec.dir = ::forward;

	// handle another boundary case...
	if (abs(asin) < ZERO_CUTOFF)
	{
		tVec.pos = sVec.x + 0.5*size.x;
		sVec.x = -0.5*size.x;
		totalLength = size.x;
		return;
	}
	

	double newX;
	double newY;
	double maxLength;
	
	newX = -0.5*size.x;
	newY = sVec.y - (sVec.x + 0.5*size.x) * (asin / acos);
	if (newY < -0.5*size.y)
	{
		newX -= (newY + 0.5*size.y) * acos / asin;
		newY = -0.5*size.y;
	}
	else if (newY > 0.5*size.y)
	{
		newX -= (newY - 0.5*size.y)* acos / asin;
		newY = 0.5*size.y;
	}
	tVec.pos = sqrt(pow(sVec.x - newX,2) + pow(sVec.y - newY,2));
	sVec.x = newX;
	sVec.y = newY;
	
	if (asin > 0) // positive slope
	{
		maxLength = (0.5*size.x - newX)/acos;
		if (newY + maxLength * asin > 0.5*size.y)
			maxLength = (0.5*size.y - newY)/asin;
	}
	else	// negative slope
	{
		maxLength = (0.5*size.x - newX)/acos;
		if (newY + maxLength * asin < -0.5*size.y)
			maxLength = (-0.5*size.y - newY)/asin;

	}
	totalLength = maxLength;

	#ifdef DBG_ACID_TEST
	if (totalLength <= 0)
	{
		cerr << "DBG/ACID TEST: calculated trajectory length is smaller than zero. [region= rectangle, x,y=" \
		<< sVec.x << ", " << sVec.y << ", angle=" << sVec.angle << ", length=" << maxLength << "]\n";
		
		
	}
	#endif

	return;
}

#ifndef NO_INLINE 
inline 
#endif
RectangleSide GRectangle::locateSide(double& x, double& y)
{
	if (abs(abs(x) - 0.5*size.x) < abs(abs(y) - 0.5*size.y))
	{
		// closer to left or right edge
		if (x > 0)
			return s_right;
		else
			return s_left;
	}
	else
	{
		// closer to top or bottom
		if (y > 0)
			return s_top;
		else
			return s_bottom;
	}
}



/*********************************** DISC REGION FUNCTIONS *******************************/

Disc::Disc(double r, Geometry* g, int idx, RegionAccountingType acc, RegionParameterType par) : 
		Cartesian(r_disc, g, idx, PI*r*r, acc, par), 
		radius(r)
{
	return;
}


SurfaceVector Disc::randomSurfaceVector()
{
	#ifdef DBG_GEOMETRY
	cout << "DBG/GEOMETRY: Disc::randomSurfaceVector() called\n";
	#endif
	SurfaceVector v;
	System* s = geometry->system;
	do
	{
		// make sure that no microtubules are nucleated exactly on a region boundary
		// this can lead to event ordering problems
		v.x = s->randomGen.randDblExc(2*radius) - radius;
		v.y = s->randomGen.randDblExc(2*radius) - radius;
	} while (pow(v.x,2) + pow(v.y,2) >= pow(radius,2));
	v.angle = s->randomGen.randExc(2*PI);
	v.region = this;

	return v;	
}

void Disc::getTrajectoryCoordinates(SurfaceVector& sVec, double& totalLength, TrajectoryVector& tVec)
{
	#ifdef DBG_ACID_TEST
	if (pow(sVec.x,2) + pow(sVec.y,2) > pow(radius + ZERO_CUTOFF,2))
	{
		cerr << "Disc::getTrajectoryCoordinates() received an out of bounds coordinate.\n";
		exit(-1);
	}
	#endif

	tVec.trajectory = NULL;
	double acos = cos(sVec.angle);
	double asin = sin(sVec.angle);

	if (acos<0)
	{
		acos *= -1;
		asin *=-1;
		sVec.angle += PI;
		tVec.dir = backward;
		if (sVec.angle>2*PI)
			sVec.angle-= 2*PI;
	}
	else
		tVec.dir = ::forward;

	// rotate -angle
	// calculate length and base
	// rotate angle
	
	double newX = acos*sVec.x + asin*sVec.y;
	double newY = -asin*sVec.x + acos*sVec.y;
	double temp = sqrt(radius*radius - newY*newY);
	
	totalLength = 2*temp;
	tVec.pos = newX + temp;
	newX = -temp;

	sVec.x = acos*newX - asin*newY;
	sVec.y = asin*newX + acos*newY;

	return;
}


