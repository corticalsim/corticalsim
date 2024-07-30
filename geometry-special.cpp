#include <fstream>
#include <iostream>
#include <iomanip>
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

void Periodic::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
        cerr <<  "Periodic::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) not implemented\n";
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

void Periodic::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  // NOT IMPLEMENTED
  return;
}

void Periodic::fixRegion(SurfaceVector & s) {
  // NOT IMPLEMENTED
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
        if (angle >= PI)        // include it twice, because you *can* start with 2PI + epsilon
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

Grid::Grid(Coord2D inSize, int n, double gridAspectratio, System* s) : 
    Geometry(g_grid, s, inSize.x*inSize.y), 
    size(inSize),
    xNumber(n)
{
    double gridsize_x = inSize.x/n;
    double gridsizeTarget_y = gridAspectratio * gridsize_x;
    yNumber= static_cast <int> (inSize.y/(gridsizeTarget_y - ZERO_CUTOFF));
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

void Grid::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
    double orientation[3][2] = {{1, 0},{0,1},{0,0}};
    int nBins,xBins,yBins,firstBin,lastBin,spareBin;
    double * cuts;
    double dx = size.x/xNumber;
    double dy = size.y/yNumber;
    double xPos, yPos,*areas;
    int i,j,k;
    sm.reset(); // set measurement to 0;
    xBins=yBins=0;
    if (sm.ori == o_x) 
    {
        xBins = sm.bins/xNumber + 4;
        cuts = new double[xBins ] ; // "certainly enough"
        areas = new double[xBins ] ; // "certainly enough"
    }
    else if (sm.ori == o_y)
    {
        yBins = sm.bins/yNumber + 4;
        cuts = new double[yBins ] ; // "certainly enough"
        areas = new double[yBins ] ; // "certainly enough"
        //cuts = new double[sm.bins/yNumber + 4 ] ; // "certainly enough"
    }
    cuts[0] = -VERY_LARGE;

    yPos = 0.5*size.y - dy;
    if (sm.ori == o_y)
        firstBin = lastBin = sm.bins-1;
    for (i=0 ; i < yNumber ; i++) {
        if (sm.ori == o_x)
            firstBin = lastBin = 0;
        xPos = -0.5*size.x + dx;
        if (sm.ori == o_y) {
            while (sm.cuts[firstBin] > yPos)
            {
                firstBin -= 1;
                if (firstBin < 0)
                {
                    break;
                }
            }
            firstBin += 1;
            nBins = lastBin-firstBin+1;
            //for (k=0 ; k < nBins ; k++)
            for (k=0 ; k < yBins-1 ; k++) {
                if ( k >= nBins)
                    cuts[k+1] = VERY_LARGE;
                else
                    cuts[k+1] = sm.cuts[k+firstBin] - yPos - 0.5*dy;
                areas[k] = dx*( min(0.5*dy, cuts[k+1]) - max(-0.5*dy,min(0.5*dy,cuts[k])) );
            }
            //cuts[k] = sm.cuts[k+firstBin] - yPos - 0.5*dy;
            lastBin = firstBin;
            if (abs(2*cuts[nBins] - dy) < ZERO_CUTOFF) {
                spareBin = 1;
            }
            else {
                spareBin = 0;
            }
        }
        for (j=0 ; j<xNumber ; j++ ) {
            if (sm.ori == o_x) {
                while (sm.cuts[lastBin] < xPos)
                {
                    lastBin += 1;
                    if (lastBin >= sm.bins )
                    {
                        break;
                    }
                }
                //lastBin -= 1;
                nBins = lastBin-firstBin+1;
                for (k=0 ; k < xBins-1 ; k++)
{
                    if ( k >= nBins)
                        cuts[k+1] = VERY_LARGE;
                    else
                        cuts[k+1] = sm.cuts[k+firstBin] -xPos + 0.5*dx;
                    //cuts[k] = sm.cuts[k+firstBin] -xPos + 0.5*dx;
                areas[k] = dy*( min(0.5*dx, cuts[k+1]) - max(-0.5*dx,min(0.5*dx,cuts[k])) );
}
                if (abs(2*cuts[nBins] - dx) < ZERO_CUTOFF) {
                    //cout << "boundary!\n";
                    spareBin = 1;
                }
                else {
                    spareBin = 0;
                }
            }
#ifdef DBG_SPATIAL
            cout << firstBin << "\t" << lastBin << "\t" << nBins << "\t" << xBins << "\t" << yBins << "\t" << xPos << "\t" << xPos+dx << "\t" << yPos << "\t" << yPos+dy << "\t" << sm.cuts[firstBin] << "\t" << sm.cuts[firstBin+nBins-1] <<  "\t" << cuts[0] << "\t" << cuts[1] << "\t" << cuts[nBins] << "\t" << cuts[nBins+1] << "\t" << OrientationTypeText[sm.ori] << "\n";
#endif
            static_cast<Cartesian*>(regions[i*xNumber+j])->getOrderParametersRawFlatBinned(sm, orientation,cuts, firstBin, nBins, spareBin, sm.ori,areas);
            xPos += dx;
            if (sm.ori == o_x )
            {
                firstBin = lastBin;
                //if (cuts[firstBin] < xPos ){
                    //firstBin +=1;
                //}
            }
        }
        yPos -= dy;
    }
#ifdef DBG_SPATIAL
    if (sm.ori == o_x) 
    {
        double nCheck = 0;
        for (k=0 ; k<sm.bins ; k++){
            nCheck += sm.nHist[k]; 
        }
        cout << "X-count: " << nCheck << "\n";
    }
    else if (sm.ori == o_y)
    {
        double nCheck = 0;
        for (k=0 ; k<sm.bins ; k++){
            nCheck += sm.nHist[k]; 
        }
        cout << "Y-count: " << nCheck << "\n";
    }
#endif
    
    sm.process(); // Important: process, writeToFile and reset must be called, in this order! (or reset before start of measurement)
    sm.writeToFile(t);

    delete [] cuts; 
    delete [] areas; 
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

void Grid::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  double xMax = 0.5* size.x;
  double yMax = 0.5* size.y;
  int regionIndex;
  regionIndex = s.region->geometryRegionIndex;
  int regionColumn = (regionIndex)%xNumber;
  int regionRow = (regionIndex)/xNumber;
  double xBase = (regionColumn-(xNumber-1)/2.)*size.x/xNumber;
  double yBase = -(regionRow-(yNumber-1)/2.)*size.y/yNumber;
  if ( s.x + xBase + x1 > xMax || s.x + xBase + x1 < -xMax ) {
    if (s.x + xBase + x2 < -xMax || s.x + xBase + x2 > xMax ) {
      cout << "Grid::shiftSurfaceVectorNoWrap: Illegal shift options X: " << x1 << ", " << x2 << ". Not shifting. \n";
    }
    else {
      s.x += x2;
    }
  } 
  else {
    s.x += x1;
  }
  if ( s.y + yBase + y1 > yMax || s.y + yBase + y1 < -yMax ) {
    if (s.y + yBase + y2 < -yMax || s.y + yBase + y2 > yMax ) {
      cout << "Grid::shiftSurfaceVectorNoWrap: Illegal shift options Y: " << y1 << ", " << y2 << ". Not shifting. \n";
    }
    else {
      s.y += y2;
    }
  } 
  else {
    s.y += y1;
  }
  fixRegion(s);
  return;
}

void Grid::fixRegion(SurfaceVector & s) {
  Region* rOld;
  double xMax, yMax; 
  rOld = s.region;
  // raise error if rOld->type != rectangle ; superfluous on Grid!
  if (rOld->type != r_rectangle ) {
    cerr << "Grid::fixRegion only implemented for rectangular parts. Exit.\n";
    exit(-32);
  }
  xMax = 0.5* size.x/xNumber;
  yMax = 0.5* size.y/yNumber;
  //yMax = (static_cast<GRectangle*>rOld) -> size.y;
  if (abs(s.x) < xMax && abs(s.y) < yMax) {
    return;
  }
  int regionIndex;
  regionIndex = rOld->geometryRegionIndex;
  // WARNING: really ugly - improve at a later stage!!
  //  for (regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);
  int regionColumn = regionIndex%xNumber;
  int regionRow = regionIndex/xNumber;
  //cout << regionRow << " " << regionColumn << " " << regionRow*xNumber + regionColumn << " " << xNumber << " " << yNumber << "\n";
  while (s.x > 0 && s.x > xMax) {
    regionColumn ++;
    s.x -= 2* xMax;
    if (regionColumn >= xNumber )
      regionColumn -= xNumber;
  }
  //cout << regionColumn<< "\n";
  while (s.x <= 0 && s.x < -xMax) {
    regionColumn --;
    s.x += 2* xMax;
    if (regionColumn <  0 )
      regionColumn += xNumber;
  }
  //cout << regionColumn<< "\n";
  if (s.x == xMax) // prevent locations exactly at the boundary
    s.x -= ZERO_CUTOFF;
  while (s.y > 0 && s.y > yMax) {
    regionRow ++;
    s.y -= 2* yMax;
    if (regionRow >= yNumber )
      regionRow -= yNumber;
  }
  while (s.y <= 0 && s.y < -yMax) {
    regionRow --;
    s.y += 2* yMax;
    if (regionRow < 0 )
      regionRow += yNumber;
  }
  if (s.y == yMax) // prevent locations exactly at the boundary
    s.y -= ZERO_CUTOFF;
  //cout << regionRow << " " << regionColumn << " " << regionRow*xNumber + regionColumn << "\n";
  s.region = regions [regionRow*xNumber + regionColumn];
  return;
}

TrajectoryVector Grid::extendTrajectory(Trajectory* oldtr, Direction dir)
{
#ifdef DBG_GEOMETRY
    cout << "DBG/GEOMETRY: Grid::extendTrajectory() called\n";
#endif
#ifdef DBG_ASTER
    cout << "DBG/ASTER: Grid::extendTrajectory() called with Direction " << dir << "\n";
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
    //  for (regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);
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
            if (angle >= PI)        // include it twice, because you *can* start with 2PI + epsilon
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

void Wormhole::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
        cerr <<  "Wormhole::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) not implemented\n";
        return;
}

double Wormhole::calculateHistogram(double hist[], bool opt)
{
    return  regions[0]->calculateHistogram(hist, opt);

}

void Wormhole::outputSnapshot(ostream& out)
{
    out << "canvas rectangle " << size.x << " " << size.y << "\n";
    out << "base 0 0 10 1 0 0 0 1 0\n";
    regions[0]->outputSnapshot(out);
    return;
}

void Wormhole::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  // NOT IMPLEMENTED
  return;
}

void Wormhole::fixRegion(SurfaceVector & s) {
  // NOT IMPLEMENTED
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

void Cylinder::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
        cerr <<  "Cylinder::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) not implemented\n";
        return;
}

double Cylinder::calculateHistogram(double hist[], bool opt)
{
    return  regions[0]->calculateHistogram(hist, opt);

}

void Cylinder::outputSnapshot(ostream& out)
{
    out << "canvas cylinder " << length << " " << radius << "\n";
    out << "base 0 0 0 1 0 0 0 1 0\n";
    regions[0]->outputSnapshot(out);
    out << "canvas disc " << radius << "\n";        // left cap (upside down!)
    out << "base " << -0.5*length << " 0 0 0 0 -1 0 -1 0\n";
    regions[1]->outputSnapshot(out);
    out << "canvas disc " << radius << "\n";        // right cap
    out << "base " << 0.5*length << " 0 0 0 0 -1 0 1 0\n";
    regions[2]->outputSnapshot(out);
    return;
}

void Cylinder::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  // NOT IMPLEMENTED
  return;
}

void Cylinder::fixRegion(SurfaceVector & s) {
  // NOT IMPLEMENTED
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
                cosEdgeAngle = 1.0; // smooth transition
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

GridCylinder::GridCylinder(double len, double r, int n, double gridAspectratio, System* s) : 
    Geometry(g_gridcylinder, s, 2*PI*r*(len + r)), 
    radius(r),
    length(len),
    number(n)
{
    gridLenSize = len/n;
    double gridRSizeTarget = gridAspectratio * gridLenSize;
    gridRNumber = static_cast <int> (2*PI*r/ gridRSizeTarget);
    gridRSize = 2*PI*r / gridRNumber;

    regions.push_back(new Disc(r, this, 0, r_accounting_dontcount, r_param_modified));  //left cap
    regions.push_back(new Disc(r, this, 1, r_accounting_dontcount, r_param_modified));  //right cap
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

// order parameters per gridpoint
void GridCylinder::getLocalOrderParameters(vector<OrderParameters> &vops)
{
    for(int i=2; i < regions.size(); i++)
    {
        OrderParameters ops;
        OrderParametersRaw opRaw;
        int regionColumn = (i-2)%number;
        int regionRow = (i-2)/number;
        static_cast<Cartesian*>(regions[i])->getOrderParametersRawCylinder(opRaw, radius, (regionColumn-(number-1)/2.)*gridLenSize, -(regionRow-(gridRNumber-1)/2.)*gridRSize);
    
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
        ops.S2 = sqrt(opRaw.si2*opRaw.si2 + opRaw.co2*opRaw.co2);
        ops.S4 = sqrt(opRaw.si4*opRaw.si4 + opRaw.co4*opRaw.co4);
        ops.S2angle = 0.5*atan2(opRaw.si2,opRaw.co2);
        ops.S4angle = 0.25*atan2(opRaw.si4,opRaw.co4);
        ops.S2Opt = sqrt(opRaw.si2Opt*opRaw.si2Opt + opRaw.co2Opt*opRaw.co2Opt);
        ops.S4Opt = sqrt(opRaw.si4Opt*opRaw.si4Opt + opRaw.co4Opt*opRaw.co4Opt);
        ops.S2angleOpt = 0.5*atan2(opRaw.si2Opt ,opRaw.co2Opt );
        ops.S4angleOpt = 0.25*atan2(opRaw.si4Opt ,opRaw.co4Opt );
        ops.localL = opRaw.localL;
        ops.localLOpt = opRaw.localLOpt;

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

        ops.R = opRaw.extractR(ops.Rdirector);
        vops.push_back(ops);
    }

    return;
}

void GridCylinder::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
    // Unroll cylinder wall for measurement; ignore caps!
    double orientation[3][2] = {{1, 0},{0,1},{0,0}};
    int nBins,xBins,yBins,firstBin,lastBin,spareBin;
    double * cuts;
    double dx = gridLenSize ; //size.x/xNumber;
    double dy = gridRSize ; //size.y/yNumber;
    double xPos, yPos,*areas;
    int i,j,k;
    sm.reset(); // set measurement to 0;
    xBins=yBins=0;
    if (sm.ori == o_x) 
    {
        xBins = sm.bins/number + 4;
        cuts = new double[xBins ] ; // "certainly enough"
        areas = new double[xBins ] ; // "certainly enough"
    }
    else if (sm.ori == o_y)
    {
        yBins = sm.bins/gridRNumber + 4;
        cuts = new double[yBins ] ; // "certainly enough"
        areas = new double[yBins ] ; // "certainly enough"
        //cuts = new double[sm.bins/yNumber + 4 ] ; // "certainly enough"
    }
    cuts[0] = -VERY_LARGE;

    //yPos = 0.5*size.y - dy;
    yPos = PI*radius - dy;
    if (sm.ori == o_y)
        firstBin = lastBin = sm.bins-1;
    //for (i=0 ; i < yNumber ; i++) {
    for (i=0 ; i < gridRNumber ; i++) {
        if (sm.ori == o_x)
            firstBin = lastBin = 0;
        //xPos = -0.5*size.x + dx;
        xPos = -0.5*length + dx;
        if (sm.ori == o_y) {
            while (sm.cuts[firstBin] > yPos)
            {
                firstBin -= 1;
                if (firstBin < 0)
                {
                    break;
                }
            }
            firstBin += 1;
            nBins = lastBin-firstBin+1;
            //for (k=0 ; k < nBins ; k++)
            for (k=0 ; k < yBins-1 ; k++) {
                if ( k >= nBins)
                    cuts[k+1] = VERY_LARGE;
                else
                    cuts[k+1] = sm.cuts[k+firstBin] - yPos - 0.5*dy;
                areas[k] = dx*( min(0.5*dy, cuts[k+1]) - max(-0.5*dy,min(0.5*dy,cuts[k])) );
            }
            //cuts[k] = sm.cuts[k+firstBin] - yPos - 0.5*dy;
            lastBin = firstBin;
            if (abs(2*cuts[nBins] - dy) < ZERO_CUTOFF) {
                spareBin = 1;
            }
            else {
                spareBin = 0;
            }
        }
        //for (j=0 ; j<xNumber ; j++ ) {
        for (j=0 ; j<number ; j++ ) {
            if (sm.ori == o_x) {
                while (sm.cuts[lastBin] < xPos)
                {
                    lastBin += 1;
                    if (lastBin >= sm.bins )
                    {
                        break;
                    }
                }
                //lastBin -= 1;
                nBins = lastBin-firstBin+1;
                for (k=0 ; k < xBins-1 ; k++)
{
                    if ( k >= nBins)
                        cuts[k+1] = VERY_LARGE;
                    else
                        cuts[k+1] = sm.cuts[k+firstBin] -xPos + 0.5*dx;
                    //cuts[k] = sm.cuts[k+firstBin] -xPos + 0.5*dx;
                areas[k] = dy*( min(0.5*dx, cuts[k+1]) - max(-0.5*dx,min(0.5*dx,cuts[k])) );
}
                if (abs(2*cuts[nBins] - dx) < ZERO_CUTOFF) {
                    //cout << "boundary!\n";
                    spareBin = 1;
                }
                else {
                    spareBin = 0;
                }
            }
#ifdef DBG_SPATIAL
            cout << firstBin << "\t" << lastBin << "\t" << nBins << "\t" << xBins << "\t" << yBins << "\t" << xPos << "\t" << xPos+dx << "\t" << yPos << "\t" << yPos+dy << "\t" << sm.cuts[firstBin] << "\t" << sm.cuts[firstBin+nBins-1] <<  "\t" << cuts[0] << "\t" << cuts[1] << "\t" << cuts[nBins] << "\t" << cuts[nBins+1] << "\t" << OrientationTypeText[sm.ori] << "\n";
#endif
            //static_cast<Cartesian*>(regions[i*xNumber+j])->getOrderParametersRawFlatBinned(sm, orientation,cuts, firstBin, nBins, spareBin, sm.ori,areas);
            static_cast<Cartesian*>(regions[i*number+j+2])->getOrderParametersRawFlatBinned(sm, orientation,cuts, firstBin, nBins, spareBin, sm.ori,areas);
            xPos += dx;
            if (sm.ori == o_x )
            {
                firstBin = lastBin;
                //if (cuts[firstBin] < xPos ){
                    //firstBin +=1;
                //}
            }
        }
        yPos -= dy;
    }
#ifdef DBG_SPATIAL
    if (sm.ori == o_x) 
    {
        double nCheck = 0;
        for (k=0 ; k<sm.bins ; k++){
            nCheck += sm.nHist[k]; 
        }
        cout << "X-count: " << nCheck << "\n";
    }
    else if (sm.ori == o_y)
    {
        double nCheck = 0;
        for (k=0 ; k<sm.bins ; k++){
            nCheck += sm.nHist[k]; 
        }
        cout << "Y-count: " << nCheck << "\n";
    }
#endif
    
    sm.process(); // Important: process, writeToFile and reset must be called, in this order! (or reset before start of measurement)
    sm.writeToFile(t);

    delete [] cuts; 
    delete [] areas; 
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
    out << "canvas disc " << radius << "\n";        // left cap (upside down!)
    out << "base " << -0.5*length << " 0 0 0 0 -1 0 -1 0\n";
    regions[0]->outputSnapshot(out);
    out << "canvas disc " << radius << "\n";        // right cap
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

void GridCylinder::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  double xMax = 0.5* length; 
  double yMax = PI*radius; 
  int regionIndex;
  regionIndex = s.region->geometryRegionIndex;
  int regionColumn = (regionIndex-2)%number;
  int regionRow = (regionIndex-2)/number;
  double xBase = (regionColumn-(number-1)/2.)*gridLenSize;
  double yBase = -(regionRow-(gridRNumber-1)/2.)*gridRSize;
  if ( s.x + xBase + x1 > xMax || s.x + xBase + x1 < -xMax ) {
    if (s.x + xBase + x2 < -xMax || s.x + xBase + x2 > xMax ) {
      cout << "GridCylinder::shiftSurfaceVectorNoWrap: Illegal shift options X: " << x1 << ", " << x2 << ". Not shifting. \n";
    }
    else {
      s.x += x2;
    }
  } 
  else {
    s.x += x1;
  }
  if ( s.y + yBase + y1 > yMax || s.y + yBase + y1 < -yMax ) {
    if (s.y + yBase + y2 < -yMax || s.y + yBase + y2 > yMax ) {
      cout << "GridCylinder::shiftSurfaceVectorNoWrap: Illegal shift options Y: " << y1 << ", " << y2 << ". Not shifting. \n";
    }
    else {
      s.y += y2;
    }
  } 
  else {
    s.y += y1;
  }
  fixRegion(s);
  return;
}

void GridCylinder::fixRegion(SurfaceVector & s) {
  Region* rOld;
  double xMax, yMax; 
  rOld = s.region;
  // raise error if rOld->type != rectangle 
  if (rOld->type != r_rectangle ) {
    cerr << "GridCylinder::fixRegion only implemented for rectangular parts. Exit.\n";
    exit(-32);
  }
  xMax = 0.5* length/number;
  yMax = PI*radius/gridRNumber;
  //yMax = (static_cast<GRectangle*>rOld) -> size.y;
  if (abs(s.x) < xMax && abs(s.y) < yMax) {
    return;
  }
  int regionIndex;
  regionIndex = rOld->geometryRegionIndex;
  // WARNING: really ugly - improve at a later stage!!
  //  for (regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);
  int regionColumn = (regionIndex-2)%number;
  int regionRow = (regionIndex-2)/number;
  //cout << regionRow << " " << regionColumn << " " << regionRow*number + regionColumn << " " << number << " " << gridRNumber << "\n";
  while (s.x > 0 && s.x > xMax) {
    regionColumn ++;
    s.x -= 2* xMax;
    if (regionColumn >= number )
      regionColumn -= number;
  }
  //cout << regionColumn<< "\n";
  while (s.x <= 0 && s.x < -xMax) {
    regionColumn --;
    s.x += 2* xMax;
    if (regionColumn <  0 )
      regionColumn += number;
  }
  //cout << regionColumn<< "\n";
  if (s.x == xMax) // prevent locations exactly at the boundary
    s.x -= ZERO_CUTOFF;
  while (s.y > 0 && s.y > yMax) {
    regionRow ++;
    s.y -= 2* yMax;
    if (regionRow >= gridRNumber )
      regionRow -= gridRNumber;
  }
  while (s.y <= 0 && s.y < -yMax) {
    regionRow --;
    s.y += 2* yMax;
    if (regionRow < 0 )
      regionRow += gridRNumber;
  }
  if (s.y == yMax) // prevent locations exactly at the boundary
    s.y -= ZERO_CUTOFF;
  //cout << regionRow << " " << regionColumn << " " << regionRow*number + regionColumn << "\n";
  s.region = regions [regionRow*number + regionColumn + 2];
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
        //      for (regionIndex=2 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);
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

void Pancake::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
        cerr <<  "Pancake::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) not implemented\n";
        return;
}

double Pancake::calculateHistogram(double hist[], bool opt)
{
    return  regions[0]->calculateHistogram(hist, opt);

}

void Pancake::outputSnapshot(ostream& out)
{
    out << "canvas disc " << radius << "\n";
    out << "base 0 0 0 1 0 0 0 1 0\n";
    regions[0]->outputSnapshot(out);
    return;
}

void Pancake::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  // NOT IMPLEMENTED
  return;
}

void Pancake::fixRegion(SurfaceVector & s) {
  // NOT IMPLEMENTED
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
    regions.push_back(new GRectangle(Coord2D(x,z), this, 0, r_accounting_dontcount, r_param_modified)); // top
    regions.push_back(new GRectangle(Coord2D(x,y), this, 1));   // front
    regions.push_back(new GRectangle(Coord2D(z,y), this, 2));   // right
    regions.push_back(new GRectangle(Coord2D(x,y), this, 3));   // back
    regions.push_back(new GRectangle(Coord2D(z,y), this, 4));   // left
    regions.push_back(new GRectangle(Coord2D(x,z), this, 5, r_accounting_dontcount, r_param_modified)); // bottom
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

void Box::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
        cerr <<  "Box::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) not implemented\n";
        return;
}


double Box::calculateHistogram(double hist[], bool opt)
{
    return  regions[1]->calculateHistogram(hist, opt);

}



void Box::outputSnapshot(ostream& out)
{
    out << "canvas rectangle " << x << " " << z << "\n";    // top
    out << "base 0 " << 0.5*y << " 0 1 0 0 0 0 -1\n";
    regions[0]->outputSnapshot(out);
    out << "canvas rectangle " << x << " " << y << "\n";    // front
    out << "base 0 0 " << 0.5*z << " 1 0 0 0 1 0\n";
    regions[1]->outputSnapshot(out);
    out << "canvas rectangle " << z << " " << y << "\n";    // right
    out << "base " << 0.5*x << " 0 0 0 0 -1 0 1 0\n";
    regions[2]->outputSnapshot(out);
    out << "canvas rectangle " << x << " " << y << "\n";    // back
    out << "base 0 0 " << -0.5*z << " -1 0 0 0 1 0\n";
    regions[3]->outputSnapshot(out);
    out << "canvas rectangle " << z << " " << y << "\n";    // left
    out << "base " << -0.5*x << " 0 0 0 0 1 0 1 0\n";
    regions[4]->outputSnapshot(out);
    out << "canvas rectangle " << x << " " << z << "\n";    // bottom
    out << "base 0 " << -0.5*y << " 0 1 0 0 0 0 1\n";
    regions[5]->outputSnapshot(out);
    return;
}

void Box::shiftSurfaceVectorNoWrap(SurfaceVector & s,double x1, double x2, double y1, double y2) {
  // NOT IMPLEMENTED
  return;
}
void Box::fixRegion(SurfaceVector & s) {
  // NOT IMPLEMENTED
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

void Cartesian::transformCuts(double *cutsIn, int binsIn, double * cutsOut, const SurfaceVector & sVec, OrientationType mOri ) 
{
    double triAngle;
    if ( mOri == o_x ) {
        triAngle = cos(sVec.angle);
        if (abs (triAngle < ZERO_CUTOFF)) {
            for (int i = 0; i<binsIn ; i++ ) {
                if (cutsIn[i] < sVec.x ) 
                    //RECENT EDIT; undone
                    cutsOut[i] = -1;            
                    //cutsOut[i] = VERY_LARGE;
                else 
                    //cutsOut[i] = -1;            
                    //RECENT EDIT; undone
                    cutsOut[i] = VERY_LARGE;
            }   
        }
        else for (int i=0 ; i<binsIn ; i++ ) {
            cutsOut[i] = (cutsIn[i] - sVec.x) / triAngle ; 
        }
    }
    if (mOri == o_y ) {
        triAngle = sin(sVec.angle);
        if (abs (triAngle < ZERO_CUTOFF)) {
            for (int i = 0; i<binsIn ; i++ ) {
                if (cutsIn[i] < sVec.y ) 
                    //RECENT EDIT ; undone
                    cutsOut[i] = -1;            
                    //cutsOut[i] = VERY_LARGE;
                else 
                    //RECENT EDIT ; undone
                    cutsOut[i] = VERY_LARGE;
                    //cutsOut[i] = -1;            
            }   
        }
        else for (int i=0 ; i<binsIn ; i++ ) {
            cutsOut[i] = (cutsIn[i] - sVec.y) / triAngle ; 
        }
    }
    //cout << "xxx\t" <<  sVec.angle << "\t" << sVec.x << "\t" << sVec.y << "\n";
    //for (int i=0 ; i< binsIn ; i++){
            //cout << cutsIn[i] << "\t" << cutsOut[i] << "\n";
    //}
    if (mOri == o_z ) { 
        cerr << "function  Cartesian::transformCuts: illegal orientation " << OrientationTypeText[mOri] << "\n";
        exit(352);
    }
    return;
}

//void Cartesian::getOrderParametersRawFlatBinned(OneDSpatialMeasurement& sm, double orientation[3][2],  double xOffset, double yOffset)
void Cartesian::getOrderParametersRawFlatBinned(OneDSpatialMeasurement& sm, double orientation[3][2],  double * localCuts, int firstBin, int nBins, int spareBin, OrientationType mOri,double * areas)
{
    int i,j,rev;

    double * sin2, *cos2, *cos4, *sin4, *Osin2, *Ocos2, *Ocos4, *Osin4;
    double angle,sin2a,cos2a,sin4a,cos4a;
    double * length, * Olength, *count;
    double * localLength, * OlocalLength;
    double * localCutsTrans;
    double u1,u2;
    double v[3];
    double * qxx, *qxy, *qxz, *qyy, *qyz, *qzz;
    Trajectory* tr;

    sin2 = new double[nBins+1+spareBin];
    cos2 = new double[nBins+1+spareBin];
    sin4 = new double[nBins+1+spareBin];
    cos4 = new double[nBins+1+spareBin];
    Osin2 = new double[nBins+1+spareBin];
    Ocos2 = new double[nBins+1+spareBin];
    Osin4 = new double[nBins+1+spareBin];
    Ocos4 = new double[nBins+1+spareBin];
    length = new double[nBins+1+spareBin]; 
    localLength = new double[nBins+1+spareBin];    
    Olength = new double[nBins+1+spareBin];    
    OlocalLength = new double[nBins+1+spareBin];   
        //RECENT EDIT; undone
    count = new double[nBins+1+spareBin]; 
    qxx = new double[nBins+1+spareBin];    
    qxy = new double[nBins+1+spareBin];    
    qxz = new double[nBins+1+spareBin];    
    qyy = new double[nBins+1+spareBin];    
    qyz = new double[nBins+1+spareBin];    
    qzz = new double[nBins+1+spareBin];    
    localCutsTrans = new double[nBins+1+spareBin];
    //for (i=0 ; i<nBins ; i++ ) {
    //qxx[i]=qxy[i]=qxz[i]=qyy[i]=qyz[i]=qzz[i]=0;
    //sin2[i]=cos2[i]=sin4[i]=cos4[i]=0;
    //Osin2[i]=Ocos2[i]=Osin4[i]=Ocos4[i]=0;
    //localLength[i] = 0.;
    //OlocalLength[i] = 0.;
    //}
    
    tr = trajectories.first();
    for (j=0; j<nBins+1+spareBin ; j++){
        sin2[j] = 0.;
        sin4[j] = 0.;
        cos2[j] = 0.;
        cos4[j] = 0.;
        localLength[j] = 0.;
        //RECENT EDIT ; undone
        count[j] = 0.;
        Osin2[j] = 0.;
        Osin4[j] = 0.;
        Ocos2[j] = 0.;
        Ocos4[j] = 0.;
        OlocalLength[j] = 0.;
        qxx[j] = 0.;
        qxy[j] = 0.;
        qxz[j] = 0.;
        qyy[j] = 0.;
        qyz[j] = 0.;
        qzz[j] = 0.;
    }

#ifdef DBG_SPATIAL
    //for (j=0 ; j<nBins ; j++){
    //cout << sm.opHist[j+firstBin].localL <<  " ";
//}
    //cout << "\n";

        double LdebugSum =0.;
        double RdebugSum =0.;
#endif
    while (tr != NULL)
    {
        //length = tr->segmentLength();
        angle = tr->base.angle;
        for (j=0; j<nBins+1+spareBin ; j++){
                length[j] = 0.;
                Olength[j] = 0.;
        }
        if (mOri == o_y && angle > PI)
               rev = -1;
        else 
                rev = 1;
        transformCuts (localCuts,nBins+1+spareBin,localCutsTrans,tr->base,mOri); // tr.base: angle, x, y. 
        //RECENT EDIT ; undone
        tr->coveredLengthSplit(length, Olength,count,localCutsTrans,nBins+spareBin,rev); // count broken.
        //tr->coveredLengthSplit(length, Olength,localCutsTrans,nBins+spareBin,rev);
#ifdef DBG_SPATIAL
        double debugSum =0.;
        double OdebugSum =0.;
        for(j = 0 ; j < nBins ; j++ ){
            debugSum += length[j];
            OdebugSum += Olength[j];
        }
        LdebugSum += debugSum;
        RdebugSum += tr->segmentLength();
        if ( abs(debugSum - tr->segmentLength()) > ZERO_CUTOFF) {
            for(j = 0 ; j < nBins ; j++ ){
                cout << j << " " << length[j] << " " << Olength[j] << " " << localCuts[j] << " " << localCutsTrans[j] << "\n";
            }
            cout << setprecision(15) << debugSum << " == "<< tr->segmentLength() << "; " << OdebugSum << " == " << tr->coveredLength() << " " << angle << "\n\n";
            cout << setprecision(6);
        }
#endif
        sin2a = sin(2*angle);
        sin4a = sin(4*angle);
        cos2a = cos(2*angle);
        cos4a = cos(4*angle);
        u1 = cos(angle);
        u2 = sin(angle);
        for (j=0 ; j<nBins+spareBin ; j++) {
            localLength[j] += length[j];
            //localCount[j] += count[j];
            OlocalLength[j] += Olength[j];
            sin2[j] += length[j]*sin2a;
            cos2[j] += length[j]*cos2a;
            sin4[j] += length[j]*sin4a;
            cos4[j] += length[j]*cos4a;
            Osin2[j] += Olength[j]*sin2a;
            Ocos2[j] += Olength[j]*cos2a;
            Osin4[j] += Olength[j]*sin4a;
            Ocos4[j] += Olength[j]*cos4a;

            for (i=0 ; i<3 ; i++)
                v[i] = orientation[i][0]*u1 + orientation[i][1]*u2;

            qxx[j] += length[j]*v[0]*v[0];
            qxy[j] += length[j]*v[0]*v[1];
            qxz[j] += length[j]*v[0]*v[2];
            qyy[j] += length[j]*v[1]*v[1];
            qyz[j] += length[j]*v[1]*v[2];
            qzz[j] += length[j]*v[2]*v[2];
        }
        tr = tr->next();

    }
    for (j=0 ; j<nBins+spareBin ; j++ ) {    
        sm.opHist[j+firstBin].si2 += sin2[j];
        sm.opHist[j+firstBin].si4 += sin4[j];
        sm.opHist[j+firstBin].co2 += cos2[j];
        sm.opHist[j+firstBin].co4 += cos4[j];
        sm.opHist[j+firstBin].localL += localLength[j];
        sm.opHist[j+firstBin].si2Opt += Osin2[j];
        sm.opHist[j+firstBin].si4Opt += Osin4[j];
        sm.opHist[j+firstBin].co2Opt += Ocos2[j];
        sm.opHist[j+firstBin].co4Opt += Ocos4[j];
        sm.opHist[j+firstBin].localLOpt += OlocalLength[j];
        //RECENT EDIT ; undone
        sm.nHist[j+firstBin] += count[j];

        sm.opHist[j+firstBin].Qxx += qxx[j];
        sm.opHist[j+firstBin].Qxy += qxy[j];
        sm.opHist[j+firstBin].Qxz += qxz[j];
        sm.opHist[j+firstBin].Qyy += qyy[j];
        sm.opHist[j+firstBin].Qyz += qyz[j];
        sm.opHist[j+firstBin].Qzz += qzz[j];
        
        for (i=0 ; i<3 ; i++)
        {
            sm.opHist[j+firstBin].isoWeights[i] += areas[j]*0.5*sqrt(pow(orientation[i][0],2) + pow(orientation[i][1],2));
        }
    }
#ifdef DBG_SPATIAL
    if (abs(LdebugSum - RdebugSum) > ZERO_CUTOFF) 
        cout << LdebugSum << " == " << RdebugSum << "\n";
#endif
    delete [] sin2 , cos2 , sin4 , cos4 , Osin2 , Ocos2 , Osin4 , Ocos4 , length ,  localLength ,   Olength ,   OlocalLength ,  qxx ,   qxy ,   qxz ,   qyy ,   qyz ,   qzz, localCutsTrans ;   

    return;
}

void Cartesian::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) 
{
        cerr <<  "Cartesian::getOneDmeasurement(OneDSpatialMeasurement& sm, double t ) not implemented\n";
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
    const double binsize=PI/geometry->system->p.angleHistogramBins;     // vanaf 7- 4: hoeken over -1/2 PI - 1/2 PI ( 0 in het midden! )
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
    SurfaceVector sVec;

    tr = trajectories.first();
    while (tr != NULL)
    {
        sVec = tr->base;
        // uncomment the following code to record trajectories in the snapshots
        //      out << "l " << graphics_trajectory << " ";
        //      out << sVec.x << " " << sVec.y << " " << sVec.angle << " " << tr->length;
        //      out << "\n";

        seg = tr->segments.begin();
        while (seg != tr->segments.end())
        {
            sVec = tr->base;
            sVec.x += xOffset;
            sVec.y += yOffset;
            translateVector(sVec, (**seg).start);
            if ((**seg).dir == backward)
            {
                sVec.angle -= PI;
                if (sVec.angle < 0)
                    sVec.angle += 2*PI;
            }
            out << "l " << graphics_mt << " ";
            out << sVec.x << " " << sVec.y << " " << sVec.angle << " " << (**seg).length();
            out << "\n";

            if ((**seg).isFirstInMT())
            {
                out << "p " << graphics_minus << " ";
                out << sVec.x << " " << sVec.y ;
                out << "\n";
            }
            if ((**seg).isLastInMT())
            {
                translateVector(sVec, (**seg).length());
                out << "p " << graphics_plus << " ";
                out << sVec.x << " " << sVec.y ;
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

    IntersectionItr thisItr;
    IntersectionItr thatItr;


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
                    thatItr = tr2->intersections.insert(pair<double,Intersection>(bp2,Intersection()));
                    thatItr->ISREF.otherTrajectory = tr1;

                    thisItr = tr1->intersections.insert(pair<double,Intersection>(bp1,Intersection()));
                    thisItr->ISREF.mirror = thatItr;
                    thisItr->ISREF.otherTrajectory = tr2;

                    // now update the reference in the other half of the intersection
                    thatItr->ISREF.mirror = thisItr;
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

Aster GRectangle::createAster(int nr, SurfaceVector sv)
{
	Aster aster;
	System* s = geometry->system;
	vector <double> thetas;
  	double firstTheta = s->randomGen.randExc(2*PI/nr);
  	thetas.push_back( firstTheta );
  	//cout << "b1=" << firstTheta*180/PI << endl;
  	for (int i=0; i<nr-1; i++)
  	{
  		/*if ( thetas.at(i) + 2*PI/nr > 2*PI )
  			thetas.push_back( thetas.at(i) + 2*PI/nr - 2*PI );
  		else
  			thetas.push_back( thetas.at(i) + 2*PI/nr );*/
  		/*if ( thetas.at(i) + 2*PI/nr > 2*PI )
  		{
  			cout << "Something went wront in selecting the angles of the branches. Abort!" << endl;
  			exit(1);
  		}	*/
  		thetas.push_back( thetas.at(i) + 2*PI/nr );
  	}
  	
  	aster.x = sv.x;
  	aster.y = sv.y;	
  	aster.angles = thetas;
  	aster.region = this;

  	return aster;	
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
        tVec.pos = sVec.y + 0.5*size.y;                            //size is public member of class Grid, type is Coord2D
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
    else    // negative slope
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

Aster Disc::createAster(int nr, SurfaceVector sv)
{
	Aster aster;
	System* s = geometry->system;
	vector <double> thetas;
  	double firstTheta = s->randomGen.randExc(2*PI);
  	
  	thetas.push_back( firstTheta );
  	for (int i=0; i<nr-1; i++)
  	{
  		if ( thetas.at(i) + 2*PI/nr > 2*PI )
  			thetas.push_back( thetas.at(i) + 2*PI/nr - 2*PI );
  		else
  			thetas.push_back( thetas.at(i) + 2*PI/nr );
  	}
  	
  	aster.x = sv.x;
  	aster.y = sv.y;	
  	aster.angles = thetas;
  	aster.region = this;
  		
  	return aster;	
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


