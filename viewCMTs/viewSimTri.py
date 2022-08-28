"""
Written by bandan chakrabortty @WageningenUR/AMOLF as part of PhD work (2013-2017)
"""
import os
import numpy
import mayavi.mlab
import utilityModule

######################
# visualization set-up
######################
# define geometry
geom = 'TOB'
viewMTarray = 'y' 
viewTriangles = 'n'
viewPrincipleAxis = 'n'
viewMTarrayOrientation = 'y'
geometryColor = [150,150,150,255]

# set input path                       
moviePath = './../inputData/' + geom + '/outputData' 

# start analysis 
if os.path.exists(moviePath):
    figSize = 10
    fig = mayavi.mlab.figure(size=(70*figSize,86*figSize),bgcolor=(1,1,1))
    print ('---------------------------------------------------------------------------------')
    print ('Shape type: ',geom)
    ########################
    # visualizing-geometry #
    ########################
    filenameIn = moviePath + '/ViewMeshDefault.off'
    with open(filenameIn, 'rb') as f:
        listVT = numpy.array([utilityModule.readrow(row, 4) for row in f])
    listVT = numpy.array([x for x in listVT if x != []])
    noV,noT,_ = numpy.array(listVT[0],int)
    G = numpy.reshape([gel.tolist() for gel in listVT[1:noV+1]],(-1, 3))
    Rvec_scale,centerSphere = utilityModule.sphereFit(G[:,0],G[:,1],G[:,2]) 
    CM = numpy.average(G,axis=0)
    Q = numpy.array(listVT[noV+1:listVT.size])
    Q = [gel.tolist() for gel in Q]
    Q = numpy.reshape(Q, (-1, 4))
    Q = numpy.delete(Q, 0, 1)
    rTipSize = 0.25*Rvec_scale
    Rvec_scale*=2.5
    triElements = [Q[i,:] for i in range(noT)]
    face = numpy.zeros(noT)
    filenameIn = moviePath + '/Draw_MT_segments.txt'
    fl = open(filenameIn, 'r')
    movieData = fl.readlines()                                                                                                                                                                                   
    if viewTriangles == 'y':
        mesh = mayavi.mlab.triangular_mesh(G[:,0],G[:,1],G[:,2],triElements,color=(0,0,0),line_width=2.,representation='wireframe')                      
    else:
        mesh = mayavi.mlab.triangular_mesh(G[:,0],G[:,1],G[:,2],triElements,color=(0, 0, 0),opacity=0)    
    mesh.mlab_source.dataset.cell_data.scalars = face 
    mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell data' 
    mesh.mlab_source.update()                                         
    mesh2 = mayavi.mlab.pipeline.set_active_attribute(mesh,cell_scalars='Cell data')                                                                       
    mySurf=mayavi.mlab.pipeline.surface(mesh2)        
    mySurf.module_manager.scalar_lut_manager.lut.table = numpy.array([geometryColor,geometryColor])
    ###################
    # visualizing MTs #
    ###################
    if (viewMTarray == 'y'):
        movieSourceMain = [ ]
        for i in range(len(movieData)):
            movieSourceMain.append(movieData[i].rstrip('\n'))
        movieSource = movieSourceMain[-1]
        shapeAtr,segInfo = movieSource.split('*')
        indvseg = segInfo.split('@')
        indvseg = indvseg[:-2]
        MT = []
        for seg in indvseg:
            segcorr = seg.split(',')
            segcorr = segcorr[:-1]
            segcorr[:] = [float(cor) for cor in segcorr]
            p1,p2 = utilityModule.split_list(segcorr)
            MT.append(p1)
            MT.append(p2)
        MT = numpy.reshape(MT, (-1, 3))            
        x1 = MT[:,0]
        y1 = MT[:,1]
        z1 = MT[:,2]
        segNum = MT.shape[0]
        connections = []
        connections.append(numpy.vstack([numpy.arange(0,segNum-1.5,2),numpy.arange(1,segNum,2)]).T)
        connections = numpy.vstack(connections)  
        src = mayavi.mlab.pipeline.scalar_scatter(x1,y1,z1)
        src.mlab_source.dataset.lines = connections
        lines= mayavi.mlab.pipeline.surface(src) # alternative: lines= mayavi.mlab.pipeline.tube(src,tube_radius=0.05,tube_sides=6)
        mayavi.mlab.pipeline.surface(lines,color=(0,1,0),line_width=8,opacity=1)
    ####################################
    # visualizing MT array orientation #
    ####################################
    if viewMTarrayOrientation == 'y':
        filenameIn = moviePath + '/measurements.txt'
        measurements = numpy.genfromtxt(filenameIn,comments=None,missing_values="nan",filling_values="0",dtype=float)
        endTimePoint = measurements.shape[0]-1
        endMeasurement = measurements[endTimePoint,:]
        Rvec = endMeasurement[-3:]
        Rvec = [a*Rvec_scale for a in Rvec]
        x0 = CM[0]
        y0 = CM[1]
        z0 = CM[2] 
        mayavi.mlab.points3d(x0-Rvec[0],y0-Rvec[1],z0-Rvec[2],color = (1.0,0.0,1.0), scale_factor = rTipSize)
        mayavi.mlab.quiver3d(x0,y0,z0,-Rvec[0],-Rvec[1],-Rvec[2],mode='arrow',scale_factor=1.0,scale_mode='vector',color=(1,0,1))
    ##############################
    # visualizing principle axes #
    ##############################
    if viewPrincipleAxis == 'y':
        filenameIn = moviePath + '/MeshInfo.txt'     
        meshStat = numpy.genfromtxt(filenameIn,delimiter=' ', usecols=0, dtype=str)         
        txtFloat = []
        for row in range(len(meshStat)):
            getfromMeshData = meshStat[row]                
            txtFloat.append(getfromMeshData.split('\t'))                                                                                                                                                                                                                                                                                       
        growth_axes = []     
        gaxCol = []        
        PA = []
        strPA =txtFloat[3][1].split(',')
        for i in range(3):
            PA.append(float(strPA[i]))
        gaxCol.append((0,0,1))
        growth_axes.append(PA/numpy.linalg.norm(PA))                                                                                                                       
        PB = []
        strPB =txtFloat[4][1].split(',')
        for i in range(3):
            PB.append(float(strPB[i]))
        gaxCol.append((1,0,0))
        growth_axes.append(PB/numpy.linalg.norm(PB))                                                                    
        PC = []
        strPC =txtFloat[5][1].split(',')
        for i in range(3):
            PC.append(float(strPC[i]))
        gaxCol.append((0,1,1))
        growth_axes.append(PC/numpy.linalg.norm(PC)) 
        pcol = 0
        for p in growth_axes[:]:
            paxwidth  = Rvec_scale
            mayavi.mlab.quiver3d([CM[0]],[CM[1]],[CM[2]],[p[0]],[p[1]],[p[2]],mode='arrow',scale_factor=paxwidth,scale_mode='none',color=gaxCol[pcol])
            mayavi.mlab.quiver3d([CM[0]],[CM[1]],[CM[2]],[-p[0]],[-p[1]],[-p[2]],mode='arrow',scale_factor=paxwidth,scale_mode='none',color=gaxCol[pcol])
            pcol+=1 
    print ('---------------------------------------------------------------------------------')
    mayavi.mlab.show()
else:
    print ('-----------------------------------------------------------------------------------------------------')
    print ('Input data :', moviePath , ' does not exists !')
    print ('-----------------------------------------------------------------------------------------------------')
                                
