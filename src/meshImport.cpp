#include "meshImport.h"

// initialize edge
void iniEdgeRecord(struct edgeRecord *head, vector<int>& v,int cid,int fid){
    for(int i=0;i<2;i++)
	head->key.push_back(v[i]);
	head->cell_id = cid;
	head->edge_index = fid;
	head->next =NULL;
}

// appending edge
void addEdgeRecord(struct edgeRecord *head, vector<int>& v,int cid,int fid) {
	edgeRecord *newedgeRecord = new edgeRecord;
	newedgeRecord->key= v;
	newedgeRecord->cell_id = cid;
	newedgeRecord->edge_index = fid;
	newedgeRecord->next = NULL;
	edgeRecord *cur = head;

	while(cur) {
		if(cur->next == NULL) {
			cur->next = newedgeRecord;
			return;
		}
		cur = cur->next;
	}
}

// delete from edge front
void deleteFirst (struct edgeRecord **head) {
    struct edgeRecord *tmp = *head;
    if (tmp == NULL) return;
    *head = tmp->next;
    free (tmp);
}

// display edge record
void display(struct edgeRecord *head)
{
	edgeRecord *list = head;
	int listCount (0);
	while(list) {
        cout << "-----------------------------"<<endl;
        cout << "Member: "<< listCount<<endl;
        cout << "-----------------------------"<<endl;

        cout << "The key ids: ";
        for(int i=0;i<2;i++)
        cout << list->key[i] <<" ";

        cout << endl;
		cout << "The cell_id: "<<list->cell_id<<endl;
		cout << "The face_index: "<<list->edge_index<<endl;
		list = list->next;
		listCount++;
	}
	cout << endl;
}

// bin-sorting
void edge_sort(edgeRecord** EList,int nodeMax)
{
    	vector<int> v1;
    	for(int i=0;i<2;i++)
    	v1.push_back(0);

    	for(int i=1;i>=0;i--)
    	{
        	edgeRecord** Bin;
        	Bin = new edgeRecord* [nodeMax];

        	for(int j = 0; j < nodeMax; j++)
        	{
           	 	Bin[j] = new edgeRecord;
            		iniEdgeRecord(Bin[j],v1,-1,-1);
        	}

        	edgeRecord* edge = *EList;
        	while(edge)
        	{
            		addEdgeRecord(Bin[edge->key[i]],edge->key,edge->cell_id,edge->edge_index);

            		edgeRecord* oldedge = edge;
            		edge = edge->next;
            		delete oldedge;
           		oldedge = NULL;
        	}

        	for(int j=0;j<nodeMax;j++)
        	deleteFirst(&Bin[j]);

        	*EList = new edgeRecord;
        	iniEdgeRecord(*EList,v1,-1,-1);


        	for(int j=0;j<nodeMax;j++)
        	{
            		edgeRecord *binedge = Bin[j];
           		while(binedge)
            		{
                		addEdgeRecord(*EList,binedge->key,binedge->cell_id,binedge->edge_index);
               			binedge=binedge->next;
            		}
         	}
            	deleteFirst(&(*EList));
            	delete[] Bin;
            	Bin = NULL;
    	}
}

// connectivity establishment
void find_connectivity(edgeRecord** FList,vector<elementList*> &eList)
{
    edgeRecord* edge = *FList;
    int ecount(0);
    while(edge!=NULL)
    {
        edgeRecord* edgenext =  edge->next;
        if(edgenext == NULL)
            break;

        if(edge->key == edgenext->key)
        {
            eList.push_back(new elementList());
            ecount = eList.size()-1;
            eList[ecount]->nodes[0]=edge->key[0];
            eList[ecount]->nodes[1]=edge->key[1];
            eList[ecount]->sharedElement[0] = edge->cell_id;
            eList[ecount]->sharedElement[1] = edgenext->cell_id;
            eList[ecount]->sharedEdge[0] = edge->edge_index;
            eList[ecount]->sharedEdge[1] = edgenext->edge_index;

            edge = edgenext->next;
        }
        else
            edge = edge->next;
    }
}

void pickup_shape(Geometry* g)
{
	// load input triangle mesh

        // on which cell to simulate ?
    	string cellIndex(g->system->p.geomParam);

        // get show mesh condition
    	int showMesh(g->system->p.showMesh);

        // get name of the input triangle mesh to be imported
    	string inFile(g->system->p.geometry);

        // scale the mesh, if asked
    	double sizeFactor(g->system->p.geometryScaleFactor);

        // global vertex list
    	vector<Vertics*> Gvertex;

        // all triangle elemet list
    	vector<Triangle3D*> element3D;

        // centroid/principle-axis of the cell
	Vector3d centroid(0.0,0.0,0.0);
   	Vector3d principleAxis(0.0,0.0,0.0);
        Vector3d axB(0.0,0.0,0.0);
        Vector3d axC(0.0,0.0,0.0);

        // list of triangle-edges to be detected
    	vector<elementList*> eList;

        int a0(0),a1(0),a2(0);
	int vmax(0),trimax(0),id(0),spf(0),cur(0),pol(0),shw(0);
	string str;

    	string filename;
   	double x(0),y(0),z(0),eNorm,eAng(0.0);

	if(inFile == "cubeReal")
    	{
		Gvertex.push_back(new Vertics(0.5*sizeFactor,0.5*sizeFactor,1.0*sizeFactor));//0
		Gvertex.push_back(new Vertics(0.5*sizeFactor,0.5*sizeFactor,0.0));//1

		Gvertex.push_back(new Vertics(0.5*sizeFactor,0.0,0.5*sizeFactor));//2
		Gvertex.push_back(new Vertics(0.5*sizeFactor,1.0*sizeFactor,0.5*sizeFactor));//3

		Gvertex.push_back(new Vertics(1.0*sizeFactor,0.5*sizeFactor,0.5*sizeFactor));//4
                Gvertex.push_back(new Vertics(0.0,0.5*sizeFactor,0.5*sizeFactor));//5

		Gvertex.push_back(new Vertics(0.0,0.0,1.0*sizeFactor));//6
		Gvertex.push_back(new Vertics(1.0*sizeFactor,0.0,1.0*sizeFactor));//7
		Gvertex.push_back(new Vertics(1.0*sizeFactor,1.0*sizeFactor,1.0*sizeFactor));//8
		Gvertex.push_back(new Vertics(0.0,1.0*sizeFactor,1.0*sizeFactor));//9

		Gvertex.push_back(new Vertics(0.0,0.0,0.0));//10
		Gvertex.push_back(new Vertics(1.0*sizeFactor,0.0,0.0));//11
		Gvertex.push_back(new Vertics(1.0*sizeFactor,1.0*sizeFactor,0.0));//12
		Gvertex.push_back(new Vertics(0.0,1.0*sizeFactor,0.0));//13

		// top-face
                /// 1st
		element3D.push_back(new Triangle3D());
		element3D[0]->vertexIds[0] = 0;
		element3D[0]->vertexIds[1] = 6;
		element3D[0]->vertexIds[2] = 7;
		element3D[0]->faceTag = 1;

                /// 2nd
		element3D.push_back(new Triangle3D());
		element3D[1]->vertexIds[0] = 0;
		element3D[1]->vertexIds[1] = 7;
		element3D[1]->vertexIds[2] = 8;
		element3D[1]->faceTag = 1;

                 /// 3rd
		element3D.push_back(new Triangle3D());
		element3D[2]->vertexIds[0] = 0;
		element3D[2]->vertexIds[1] = 8;
		element3D[2]->vertexIds[2] = 9;
		element3D[2]->faceTag = 1;

		/// 4th
		element3D.push_back(new Triangle3D());
		element3D[3]->vertexIds[0] = 0;
		element3D[3]->vertexIds[1] = 9;
		element3D[3]->vertexIds[2] = 6;
		element3D[3]->faceTag = 1;

		// bottom-face
                /// 1st
		element3D.push_back(new Triangle3D());
		element3D[4]->vertexIds[0] = 1;
		element3D[4]->vertexIds[1] = 10;
		element3D[4]->vertexIds[2] = 11;
		element3D[4]->faceTag = 1;

                /// 2nd
		element3D.push_back(new Triangle3D());
		element3D[5]->vertexIds[0] = 1;
		element3D[5]->vertexIds[1] = 11;
		element3D[5]->vertexIds[2] = 12;
		element3D[5]->faceTag = 1;

                 /// 3rd
		element3D.push_back(new Triangle3D());
		element3D[6]->vertexIds[0] = 1;
		element3D[6]->vertexIds[1] = 12;
		element3D[6]->vertexIds[2] = 13;
		element3D[6]->faceTag = 1;

		/// 4th
		element3D.push_back(new Triangle3D());
		element3D[7]->vertexIds[0] = 1;
		element3D[7]->vertexIds[1] = 13;
		element3D[7]->vertexIds[2] = 10;
		element3D[7]->faceTag = 1;

		// left-face
                /// 1st
		element3D.push_back(new Triangle3D());
		element3D[8]->vertexIds[0] = 2;
		element3D[8]->vertexIds[1] = 6;
		element3D[8]->vertexIds[2] = 7;
		element3D[8]->faceTag = 3;

                /// 2nd
		element3D.push_back(new Triangle3D());
		element3D[9]->vertexIds[0] = 2;
		element3D[9]->vertexIds[1] = 7;
		element3D[9]->vertexIds[2] = 11;
		element3D[9]->faceTag = 3;

                 /// 3rd
		element3D.push_back(new Triangle3D());
		element3D[10]->vertexIds[0] = 2;
		element3D[10]->vertexIds[1] = 11;
		element3D[10]->vertexIds[2] = 10;
		element3D[10]->faceTag = 3;

		/// 4th
		element3D.push_back(new Triangle3D());
		element3D[11]->vertexIds[0] = 2;
		element3D[11]->vertexIds[1] = 10;
		element3D[11]->vertexIds[2] = 6;
		element3D[11]->faceTag = 3;

		// right-face
                /// 1st
		element3D.push_back(new Triangle3D());
		element3D[12]->vertexIds[0] = 3;
		element3D[12]->vertexIds[1] = 8;
		element3D[12]->vertexIds[2] = 9;
		element3D[12]->faceTag = 3;

                /// 2nd
		element3D.push_back(new Triangle3D());
		element3D[13]->vertexIds[0] = 3;
		element3D[13]->vertexIds[1] = 9;
		element3D[13]->vertexIds[2] = 13;
		element3D[13]->faceTag = 3;

                 /// 3rd
		element3D.push_back(new Triangle3D());
		element3D[14]->vertexIds[0] = 3;
		element3D[14]->vertexIds[1] = 13;
		element3D[14]->vertexIds[2] = 12;
		element3D[14]->faceTag = 3;

		/// 4th
		element3D.push_back(new Triangle3D());
		element3D[15]->vertexIds[0] = 3;
		element3D[15]->vertexIds[1] = 12;
		element3D[15]->vertexIds[2] = 8;
		element3D[15]->faceTag = 3;

		// front-face
                /// 1st
		element3D.push_back(new Triangle3D());
		element3D[16]->vertexIds[0] = 4;
		element3D[16]->vertexIds[1] = 7;
		element3D[16]->vertexIds[2] = 8;
		element3D[16]->faceTag = 4;

                /// 2nd
		element3D.push_back(new Triangle3D());
		element3D[17]->vertexIds[0] = 4;
		element3D[17]->vertexIds[1] = 8;
		element3D[17]->vertexIds[2] = 12;
		element3D[17]->faceTag = 4;

                 /// 3rd
		element3D.push_back(new Triangle3D());
		element3D[18]->vertexIds[0] = 4;
		element3D[18]->vertexIds[1] = 12;
		element3D[18]->vertexIds[2] = 11;
		element3D[18]->faceTag = 4;

		/// 4th
		element3D.push_back(new Triangle3D());
		element3D[19]->vertexIds[0] = 4;
		element3D[19]->vertexIds[1] = 11;
		element3D[19]->vertexIds[2] = 7;
		element3D[19]->faceTag = 4;

		// back-face
                /// 1st
		element3D.push_back(new Triangle3D());
		element3D[20]->vertexIds[0] = 5;
		element3D[20]->vertexIds[1] = 6;
		element3D[20]->vertexIds[2] = 9;
		element3D[20]->faceTag = 4;

                /// 2nd
		element3D.push_back(new Triangle3D());
		element3D[21]->vertexIds[0] = 5;
		element3D[21]->vertexIds[1] = 9;
		element3D[21]->vertexIds[2] = 13;
		element3D[21]->faceTag = 4;

                 /// 3rd
		element3D.push_back(new Triangle3D());
		element3D[22]->vertexIds[0] = 5;
		element3D[22]->vertexIds[1] = 13;
		element3D[22]->vertexIds[2] = 10;
		element3D[22]->faceTag = 4;

		/// 4th
		element3D.push_back(new Triangle3D());
		element3D[23]->vertexIds[0] = 5;
		element3D[23]->vertexIds[1] = 10;
		element3D[23]->vertexIds[2] = 6;
		element3D[23]->faceTag = 4;

		for(int i=0;i<element3D.size();i++)
		element3D[i]->elementId = i;
    	}

	else if(inFile == "2D-plane")
	{
		Gvertex.push_back(new Vertics(0.0,0.0,0.0));
		Gvertex.push_back(new Vertics(0.0,1.0*sizeFactor,0.0));
		Gvertex.push_back(new Vertics(1.0*sizeFactor,1.0*sizeFactor,0.0));
		Gvertex.push_back(new Vertics(1.0*sizeFactor,0.0,0.0));
		Gvertex.push_back(new Vertics(0.5*sizeFactor,0.5*sizeFactor,0.0));
		Gvertex.push_back(new Vertics(0.5*sizeFactor,0.5*sizeFactor,-1.0));

		element3D.push_back(new Triangle3D());
		element3D[0]->vertexIds[0] = 0;
	 	element3D[0]->vertexIds[1] = 1;
		element3D[0]->vertexIds[2] = 4;

		element3D.push_back(new Triangle3D());
		element3D[1]->vertexIds[0] = 1;
		element3D[1]->vertexIds[1] = 2;
		element3D[1]->vertexIds[2] = 4;

		element3D.push_back(new Triangle3D());
		element3D[2]->vertexIds[0] = 2;
		element3D[2]->vertexIds[1] = 3;
		element3D[2]->vertexIds[2] = 4;

		element3D.push_back(new Triangle3D());
		element3D[3]->vertexIds[0] = 3;
		element3D[3]->vertexIds[1] = 0;
		element3D[3]->vertexIds[2] = 4;

		element3D.push_back(new Triangle3D());
		element3D[4]->vertexIds[0] = 0;
		element3D[4]->vertexIds[1] = 1;
		element3D[4]->vertexIds[2] = 5;

		element3D.push_back(new Triangle3D());
		element3D[5]->vertexIds[0] = 1;
		element3D[5]->vertexIds[1] = 2;
		element3D[5]->vertexIds[2] = 5;

		element3D.push_back(new Triangle3D());
		element3D[6]->vertexIds[0] = 2;
		element3D[6]->vertexIds[1] = 3;
		element3D[6]->vertexIds[2] = 5;

		element3D.push_back(new Triangle3D());
		element3D[7]->vertexIds[0] = 3;
		element3D[7]->vertexIds[1] = 0;
		element3D[7]->vertexIds[2] = 5;
	}

        else
	{
		// get input mest file path
		string inputMeshFile = g->system->p.inputDir+"/" + g->system->p.geometry + "/" + inFile+".dat";

		// open input mesh file
		fstream fp;
	    	fp.open(inputMeshFile.c_str(),ios::in);

	    	fp >> vmax >> trimax;

		// read the global vertex coordinates
	   	for(int i=0;i<vmax;i++)
	    	{
			fp >> x >> y >> z;
			Gvertex.push_back(new Vertics(0.0,0.0,0.0));

		        // scale the mesh, if asked
			Gvertex[i]->x = x*sizeFactor;
			Gvertex[i]->y = y*sizeFactor;
			Gvertex[i]->z = z*sizeFactor;
	    	}

		// read triangle informatio:
		// (a) vertex coordinates
		// (b) spcial face tag
		// (c) curvature tag
		// (d) Pcat norm factor
		// (d) average edge angle (of the curvature domain it belongs)
	    	for(int i=0;i<trimax;i++)
	    	{
			fp >> a0 >> a1 >> a2 >> spf >> cur  >> eNorm >> eAng;
			element3D.push_back(new Triangle3D());

				element3D[i]->elementId = i;
		                element3D[i]->faceTag = spf;
		                element3D[i]->edgeTag = cur;
				element3D[i]->edgAngNorm = eNorm;
				element3D[i]->pcatEdgeWeight = eAng;

		                // vertices need zero-based numbering
				element3D[i]->vertexIds[0] = a0-1;
				element3D[i]->vertexIds[1] = a1-1;
				element3D[i]->vertexIds[2] = a2-1;
	    	}

	    	fp >> str;

		// close the input mesh file
	    	fp.close();
	}

        // establish connectivity among the triangles
   	makeGraph(Gvertex.size(),eList,element3D);

        // store the numer of triangles, vertices and edges
        int F = element3D.size();
    	int V = Gvertex.size();
    	int E = eList.size();

        // check mesh consistency
        bool convex(checkEular(F,V,E));

	// copy triangle elements for viewing mesh (visualization in meshLab)
	vector<Triangle3D*>viewTriangle;
	for(int i = 0; i < element3D.size(); i++)
	{
		viewTriangle.push_back(new Triangle3D());
		viewTriangle[i]->vertexIds[0] = element3D[i]->vertexIds[0];
		viewTriangle[i]->vertexIds[1] = element3D[i]->vertexIds[1];
		viewTriangle[i]->vertexIds[2] = element3D[i]->vertexIds[2];

		viewTriangle[i]->faceTag = element3D[i]->faceTag;
		viewTriangle[i]->edgeTag = element3D[i]->edgeTag;
	}

	// check surface orientation
    	orientSurface(Gvertex,element3D,viewTriangle);

        // check surface orientation
        bool surfaceOriented(checkSurfaceOrientation(element3D));

        // measure rigidbody properties
	rigidBodyProperties(Gvertex,element3D,centroid,principleAxis,axB,axC);

        // assign coordinate values to the triangle vertices
    	connectWithGlobe(Gvertex,element3D);

        // get triangle edge related measures/properties
    	edgeDescriptors(Gvertex,eList,element3D);

	/// PBC for 2D-plane
    	if(inFile =="2D-plane")
	{
		principleAxis << 1.0,0.0,0.0;
		axB << 0.0,1.0,0.0;
		axC << 0.0,0.0,0.0;
    		establishPBC(Gvertex,element3D,viewTriangle,eList,sizeFactor);
	}

	// visualize mesh in meshLab
	viewGraph(Gvertex,viewTriangle,centroid,g->system->p.outputDir + "/" + g->system->p.geometry + "/outputData");

	double totalArea(0.0);
        for(int eno=0;eno<element3D.size();eno++)
	totalArea+=element3D[eno]->area;


        // show mesh in enabled
    	if(showMesh==1)
    	{
         	cout << "============================================"<<endl;

                // show principle axis
    		cout << "Principle Axis ---> (" << principleAxis(0,0)<<"," <<principleAxis(1,0)<<","<<principleAxis(2,0) <<")"<<endl;

                // show the numer of triangles, vertices and edges
    		cout << "F = "<< F << " V = "<< V << " E = "<< E <<endl;

                // mesh is non-convex
     		if(!convex)
     		cout << "The surface mesh is non-Euler" <<endl;

                // mesh in convex
        	else
        	cout << "The surface mesh is Euler" <<endl;

                if(surfaceOriented)
    		cout << "Surface oriented (outwards) successfully"<<endl;

		cout << "============================================"<<endl;
	}

        // summary of the mesh properties
	fstream fp;
        filename = g->system->p.outputDir + "/" + g->system->p.geometry + "/outputData" + string("/MeshInfo.txt");
	fp.open(filename.c_str(),ios::out);

	fp << "CellIndex\t" << g->system->p.geomParam <<endl;
	fp << "TriangleNumbers\t" << element3D.size() <<endl;
	fp << "Centroid\t" << centroid(0,0)<<","<< centroid(1,0)<<","<< centroid(2,0)<<endl;
	fp << "PA\t" << principleAxis(0,0) <<","<< principleAxis(1,0)<< "," << principleAxis(2,0) <<endl;
        fp << "PB\t" << axB(0,0) <<","<< axB(1,0)<< "," << axB(2,0) <<endl;
	fp << "PC\t" << axC(0,0) <<","<< axC(1,0)<< "," << axC(2,0) <<endl;
	fp << "Area\t" << totalArea <<endl;

	fp.close();

        // starting simulation on the surface-mesh
	map<string, double>::iterator it;

	// list of edge catastrophe values
        int edgNumber = g->system->p.edgCatMap.size();
	double pCatList[edgNumber];
        for(int i=0;i<= edgNumber;i++)
        pCatList[i]= 0.0;

	for(it = g->system->p.edgCatMap.begin(); it!= g->system->p.edgCatMap.end(); it++)
	{
                vector<string> currEdg = split(it->first, '_');
                int edgIdx = boost::lexical_cast<int>(currEdg[1]);
		double edgCatVal = it->second;
                pCatList[edgIdx] = edgCatVal;
	}

	// list of face catastrophe values
        int faceNumber = g->system->p.faceCatMap.size();
        for(int i=0;i<= faceNumber;i++)
        g->system->p.RegionKcatMultiplier.push_back(1.0);

	for(it = g->system->p.faceCatMap.begin(); it!= g->system->p.faceCatMap.end(); it++)
	{
                vector<string> currFace = split(it->first, '_');
                int faceIdx = boost::lexical_cast<int>(currFace[1]);
		double faceCatVal = it->second;
                g->system->p.RegionKcatMultiplier[faceIdx] = faceCatVal;
	}

        // assign the centroid and principle axis to the geometry
	g->objectCM = centroid;
	g->objectPA = principleAxis;
        g->nucleousPosition << 0.0,0.0,0.0;

        // total surface area of the geometry
	g->area = totalArea;

        // create only a single path area (at first place)
        g->patchArea[0]= g->area;

	// include all triangles to a single path (at first place)
	for(int eno=0;eno<element3D.size();eno++)
	g->RegionsIndex[0].push_back(eno);

        // total number of region in the geometry
	g->elementMax = element3D.size();

        // to create PPB
        // (a)  global vertex to be used to create PPB
        for(int i=0;i< Gvertex.size();i++)
	{
                g->globalVertex.push_back(Vertics(0.0,0.0,0.0));
		g->globalVertex[i].x = Gvertex[i]->x;
		g->globalVertex[i].y = Gvertex[i]->y;
		g->globalVertex[i].z = Gvertex[i]->z;
	}

        // (b) element list
        for(int i=0;i< eList.size();i++)
	{
		g->ppbEdgeList.push_back(elementList());

		for(int j=0;j<2;j++)
		{
    			g->ppbEdgeList[i].nodes[j] = eList[i]->nodes[j];
			g->ppbEdgeList[i].sharedElement[j] = eList[i]->sharedElement[j];
			g->ppbEdgeList[i].sharedEdge[j] = eList[i]->sharedEdge[j];
		}
	}

        // assign memory for all the 2d-triangular regions of the geometry
	for(int eno = 0; eno < element3D.size(); eno++)
	g->regions.push_back(new Triangle(element3D[eno]->area,g));

	// copy 3d-triangle properties information to 2d-triangular regions
	Image3dTo2D(g->regions,element3D,pCatList);


    	// 3d-triangle properties are transfered to 2d-triangular regions, so release all the memory occupied by 3d triangles
    	for(int eno  = 0; eno < element3D.size(); eno++)
    	{
        	delete element3D[eno];
    		element3D[eno] == NULL;
    	}
    	element3D.clear();

	// global vertices no longer needed, so release the memories occupied by the global vertices
    	for(int gno = 0; gno < Gvertex.size(); gno++)
    	{
        	delete Gvertex[gno];
    		Gvertex[gno] = NULL;
    	}
    	Gvertex.clear();

        // edge list is no longer needed, so release the memories occupied by the edges
	for(int el = 0; el < eList.size(); el++)
    	{
        	delete eList[el];
    		eList[el] = NULL;
    	}
	eList.clear();

    return;
}

void makeGraph(int Gvertex,vector<elementList*>& eList, vector<Triangle3D*>& element3D)
{
	// make the connectivity graph
	void edge_sort(edgeRecord**,int);
	void find_connectivity(edgeRecord**,vector<elementList*> &);

    	struct edgeRecord *EdgList = new edgeRecord;

    	vector<int> v;
    	v.push_back(0);
    	v.push_back(0);

    	iniEdgeRecord(EdgList,v,-1,-1);

    	v.clear();

    	for(int i=0;i<element3D.size();i++)
    	{
        	for(int j=0;j<3;j++)
        	element3D[i]->vertexIdsMap[element3D[i]->vertexIds[j]] = j;

        	for(int j=0;j<3;j++)
        	{
            		v.push_back(element3D[i]->vertexIds[j%3]);
            		v.push_back(element3D[i]->vertexIds[(j+1)%3]);
            		element3D[i]->side[j].orientation=v;

            		element3D[i]->side[j].orientationMap[v[0]]=j%3;
            		element3D[i]->side[j].orientationMap[v[1]]=(j+1)%3;
            		element3D[i]->side[j].dir.push_back(j%3);
            		element3D[i]->side[j].dir.push_back((j+1)%3);

            		element3D[i]->side[j].excludePoint = element3D[i]->vertexIds[(j+2)%3];
            		element3D[i]->side[j].excludePointMap[element3D[i]->vertexIds[(j+2)%3]] = (j+2)%3;

            		sort(v.begin(),v.end());

            		addEdgeRecord(EdgList,v,i,j);
            		v.clear();
        	}
        	element3D[i]->elementId = i;
    	}
    	deleteFirst(&(EdgList));

    	edge_sort(&EdgList,Gvertex);
    	find_connectivity(&EdgList,eList);

    	for(int i=0;i<eList.size();i++)
    	{
       	 	element3D[eList[i]->sharedElement[0]]->sideMap[eList[i]->sharedEdge[0]]=eList[i]->sharedElement[1];
        	element3D[eList[i]->sharedElement[1]]->sideMap[eList[i]->sharedEdge[1]]=eList[i]->sharedElement[0];

        	element3D[eList[i]->sharedElement[0]]->sideRevMap[eList[i]->sharedElement[1]] = eList[i]->sharedEdge[0];
        	element3D[eList[i]->sharedElement[1]]->sideRevMap[eList[i]->sharedElement[0]] = eList[i]->sharedEdge[1];

        	element3D[eList[i]->sharedElement[0]]->sideToEdge[eList[i]->sharedEdge[0]] = i;
        	element3D[eList[i]->sharedElement[1]]->sideToEdge[eList[i]->sharedEdge[1]] = i;
    	}

    	delete EdgList;
    	EdgList = NULL;

    	return;
}

void orientSurface(vector<Vertics*>& Gvertex,vector<Triangle3D*>& element3D,vector<Triangle3D*>& viewTriangle)
{
	// orient surface normals
    	void correctOrientation(vector<Vertics*>&,Triangle3D*&,Triangle3D*&,int,int,int &);

    	vector<Triangle3D*>TriangleAlive;
    	TriangleAlive.push_back(element3D[0]);

    	vector<Triangle3D*> TriangleDead;

    	while(!TriangleAlive.empty())
    	{
        	Triangle3D* nextTriangle = TriangleAlive.back();
        	TriangleAlive.pop_back();

        	TriangleDead.push_back(nextTriangle);

        	for(int i=0;i<3;i++)
        	{
            		Triangle3D* neighTriangle = element3D[nextTriangle->sideMap[i]];

            		bool isPresent=find(TriangleDead.begin(),TriangleDead.end(),neighTriangle)!=TriangleDead.end();
            		if(!isPresent)
            		{
				int flip(1);
                		int n1 = neighTriangle->sideRevMap[nextTriangle->elementId];
                		correctOrientation(Gvertex,nextTriangle,neighTriangle,i,n1,flip);
                		TriangleAlive.push_back(neighTriangle);
				if(flip == -1)
				{
					viewTriangle[nextTriangle->sideMap[i]]->vertexIds[1] = element3D[nextTriangle->sideMap[i]]->vertexIds[2];
					viewTriangle[nextTriangle->sideMap[i]]->vertexIds[2] = element3D[nextTriangle->sideMap[i]]->vertexIds[1];
				}
            		}
        	}
    	}
}

void correctOrientation(vector<Vertics*>& Gvertex,Triangle3D*& nextTriangle,Triangle3D*& neighTriangle,int n1,int n2, int &flip)
{
   if(nextTriangle->side[n1].orientation==neighTriangle->side[n2].orientation)
   {
        for(int j=0;j<3;j++)
        {
            reverse(neighTriangle->side[j].orientation.begin(), neighTriangle->side[j].orientation.end());
            reverse(neighTriangle->side[j].dir.begin(), neighTriangle->side[j].dir.end());
	    flip*=-1;
        }
   }

    return;
}

bool checkSurfaceOrientation(vector<Triangle3D*>& element3D)
{
    bool match(false);
    for(int i=0;i< element3D.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            match = false;
            int n = element3D[i]->sideMap[j];
            int s = element3D[n]->sideRevMap[i];
            vector<int> v;
            v.push_back(element3D[n]->side[s].orientation[1]);
            v.push_back(element3D[n]->side[s].orientation[0]);

            if(element3D[i]->side[j].orientation == v)
            match = true;

            if(!match)
            {
                cout <<"Surface not oriented ... emergency exit !!!"<<endl;
                exit(1);
            }
        }
    }

    return(match);
}

bool checkEular(int F,int V,int E)
{
    bool graph(((3*F/2)-E)==0);

    if(!graph)
    {
        cout <<"Inconsistency detected between the nuber of TRIANGLES and EDGES  ... emergency exit !!!"<<endl;
        exit(-1);
    }

    graph = ((F+V-E)==2);
    if(!graph)
    {
        cout <<"Euler condition violated ... emergency exit !!!"<<endl;
        exit(-1);
    }

    return(graph);
}

void rigidBodyProperties(vector<Vertics*>& Gvertex,vector<Triangle3D*>& element3D,Vector3d & cm,Vector3d & pax,Vector3d &pbx,Vector3d &pcx)
{
        // rigidBodyProperties
	vector<Vector3d> xRow1,yRow1,zRow1,xRow2,yRow2,zRow2,xRow3,yRow3,zRow3;
    	double x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3;
    	double x_2,y_2,z_2,xy,xz,yz;
    	double x_3,y_3,z_3;
    	double x_2y,x_2z,y_2x,y_2z,z_2y,z_2x,xyz;

    	for(int eno = 0; eno < element3D.size();eno++)
    	{
		// calculate normals to each triangle elements
        	x =  Gvertex[element3D[eno]->side[0].orientation[1]]->x- Gvertex[element3D[eno]->side[0].orientation[0]]->x;
        	y =  Gvertex[element3D[eno]->side[0].orientation[1]]->y- Gvertex[element3D[eno]->side[0].orientation[0]]->y;
        	z =  Gvertex[element3D[eno]->side[0].orientation[1]]->z- Gvertex[element3D[eno]->side[0].orientation[0]]->z;
        	Vector3d RL(x,y,z);

                x =  Gvertex[element3D[eno]->side[0].orientation[0]]->x- Gvertex[element3D[eno]->side[0].excludePoint]->x;
        	y =  Gvertex[element3D[eno]->side[0].orientation[0]]->y- Gvertex[element3D[eno]->side[0].excludePoint]->y;
        	z =  Gvertex[element3D[eno]->side[0].orientation[0]]->z- Gvertex[element3D[eno]->side[0].excludePoint]->z;

        	Vector3d RR(x,y,z);

        	element3D[eno]->givenNormal=RR.cross(RL);;

        	double normalMagnitude = element3D[eno]->givenNormal.norm();
        	element3D[eno]->area = 0.5*normalMagnitude;
        	element3D[eno]->givenNormal/=normalMagnitude;

		// access to element vertex
		x1 = Gvertex[element3D[eno]->vertexIds[0]]->x;
		y1 = Gvertex[element3D[eno]->vertexIds[0]]->y;
		z1 = Gvertex[element3D[eno]->vertexIds[0]]->z;

		x2 = Gvertex[element3D[eno]->vertexIds[1]]->x;
		y2 = Gvertex[element3D[eno]->vertexIds[1]]->y;
		z2 = Gvertex[element3D[eno]->vertexIds[1]]->z;

		x3 = Gvertex[element3D[eno]->vertexIds[2]]->x;
		y3 = Gvertex[element3D[eno]->vertexIds[2]]->y;
		z3 = Gvertex[element3D[eno]->vertexIds[2]]->z;

		// calculate element mid-point
		element3D[eno]->midPoint(0,0) = (x1+x2+x3)/3.0;
        	element3D[eno]->midPoint(1,0) = (y1+y2+y3)/3.0;
        	element3D[eno]->midPoint(2,0) = (z1+z2+z3)/3.0;

		// inertia utilities
        	x_2 = ((x1+x2)*(x2+x3) + x1*x1 + x3*x3)/12.0;
        	y_2 = ((y1+y2)*(y2+y3) + y1*y1 + y3*y3)/12.0;
        	z_2 = ((z1+z2)*(z2+z3) + z1*z1 + z3*z3)/12.0;

        	xy = ((x1+x2+x3)*(y1+y2+y3) + x1*y1 + x2*y2 + x3*y3)/24.0;
        	xz = ((x1+x2+x3)*(z1+z2+z3) + x1*z1 + x2*z2 + x3*z3)/24.0;
        	yz = ((y1+y2+y3)*(z1+z2+z3) + y1*z1 + y2*z2 + y3*z3)/24.0;

        	x_3=((x1+x2+x3)*(x1*x1+x2*x2+x3*x3) + x1*x2*x3)/20.0;
        	y_3=((y1+y2+y3)*(y1*y1+y2*y2+y3*y3) + y1*y2*y3)/20.0;
        	z_3=((z1+z2+z3)*(z1*z1+z2*z2+z3*z3) + z1*z2*z3)/20.0;

        	x_2y=((3.0*y1+y2+y3)*x1*x1 + (y1+3.0*y2+y3)*x2*x2 + (y1+y2+3.0*y3)*x3*x3 + (2.0*y1+2.0*y2+y3)*x1*x2 + (2.0*y1+y2+2.0*y3)*x1*x3 + (y1+2.0*y2+2.0*y3)*x2*x3)/60.0;
        	x_2z=((3.0*z1+z2+z3)*x1*x1 + (z1+3.0*z2+z3)*x2*x2 + (z1+z2+3.0*z3)*x3*x3 + (2.0*z1+2.0*z2+z3)*x1*x2 + (2.0*z1+z2+2.0*z3)*x1*x3 + (z1+2.0*z2+2.0*z3)*x2*x3)/60.0;

        	y_2x=((3.0*x1+x2+x3)*y1*y1 + (x1+3.0*x2+x3)*y2*y2 + (x1+x2+3.0*x3)*y3*y3 + (2.0*x1+2.0*x2+x3)*y1*y2 + (2.0*x1+x2+2.0*x3)*y1*y3 + (x1+2.0*x2+2.0*x3)*y2*y3)/60.0;
        	y_2z=((3.0*z1+z2+z3)*y1*y1 + (z1+3.0*z2+z3)*y2*y2 + (z1+z2+3.0*z3)*y3*y3 + (2.0*z1+2.0*z2+z3)*y1*y2 + (2.0*z1+z2+2.0*z3)*y1*y3 + (z1+2.0*z2+2.0*z3)*y2*y3)/60.0;

        	z_2y=((3.0*y1+y2+y3)*z1*z1 + (y1+3.0*y2+y3)*z2*z2 + (y1+y2+3.0*y3)*z3*z3 + (2.0*y1+2.0*y2+y3)*z1*z2 + (2.0*y1+y2+2.0*y3)*z1*z3 + (y1+2.0*y2+2.0*y3)*z2*z3)/60.0;
        	z_2x=((3.0*x1+x2+x3)*z1*z1 + (x1+3.0*x2+x3)*z2*z2 + (x1+x2+3.0*x3)*z3*z3 + (2.0*x1+2.0*x2+x3)*z1*z2 + (2.0*x1+x2+2.0*x3)*z1*z3 + (x1+2.0*x2+2.0*x3)*z2*z3)/60.0;

        	xyz=((x1+x2+x3)*(y1+y2+y3)*(z1+z2+z3) - (y2*z3+y3*z2-4.0*y1*z1)*x1/2.0 -(y1*z3+y3*z1-4.0*y2*z2)*x2/2.0 - (y1*z2+y2*z1-4.0*y3*z3)*x3/2.0)/60.0;

        	xRow1.push_back(Vector3d(x_2,2*xy,2*xz));
        	yRow1.push_back(Vector3d(2*xy,y_2,2*yz));
        	zRow1.push_back(Vector3d(2*xz,2*yz,z_2));

        	xRow2.push_back(Vector3d(x_2y,y_2x,2*xyz));
        	yRow2.push_back(Vector3d(x_2z,2*xyz,z_2x));
        	zRow2.push_back(Vector3d(2*xyz,y_2z,z_2y));

        	xRow3.push_back(Vector3d(x_3,3*x_2y,3*x_2z));
        	yRow3.push_back(Vector3d(3*y_2x,y_3,3*y_2z));
        	zRow3.push_back(Vector3d(3*z_2x,3*z_2y,z_3));
    	}

    	double m000(0.0),m100(0.0),m010(0.0),m001(0.0);
    	double m110(0.0),m101(0.0),m011(0.0),m200(0.0),m020(0.0),m002(0.0);
    	double totalArea(0.0);

    	for(int eno = 0; eno < element3D.size();eno++)
    	{
		totalArea += element3D[eno]->area;

        	// zeroth order moment (same as volume enclosed by the mesh)
        	m000 += element3D[eno]->area*element3D[eno]->givenNormal.dot(element3D[eno]->midPoint);

	        // first order moments (together they specify the centroid of the mesh)
	        m100+= element3D[eno]->area*element3D[eno]->givenNormal.dot(xRow1[eno]);
	        m010+= element3D[eno]->area*element3D[eno]->givenNormal.dot(yRow1[eno]);
	        m001+= element3D[eno]->area*element3D[eno]->givenNormal.dot(zRow1[eno]);

	        m110+= element3D[eno]->area*element3D[eno]->givenNormal.dot(xRow2[eno]);
	        m101+= element3D[eno]->area*element3D[eno]->givenNormal.dot(yRow2[eno]);
	        m011+= element3D[eno]->area*element3D[eno]->givenNormal.dot(zRow2[eno]);

	        m200+= element3D[eno]->area*element3D[eno]->givenNormal.dot(xRow3[eno]);
	        m020+= element3D[eno]->area*element3D[eno]->givenNormal.dot(yRow3[eno]);
	        m002+= element3D[eno]->area*element3D[eno]->givenNormal.dot(zRow3[eno]);
	}

	m000/=3.0;

	m100/=3.0;
    	m010/=3.0;
    	m001/=3.0;

    	m110/=3.0;
    	m101/=3.0;
    	m011/=3.0;

    	m200/=4.5;
    	m020/=4.5;
    	m002/=4.5;

    	// calculate the inertia tensor
	double Ixx,Iyy,Izz,Ixy,Ixz,Iyz;
	Ixx=(m020+m002)-((m010*m010+m001*m001)/m000);
    	Iyy=(m200+m002)-((m100*m100+m001*m001)/m000);
    	Izz=(m200+m020)-((m100*m100+m010*m010)/m000);
    	Ixy=m110-m100*m010/m000;
    	Ixz=m101-m100*m001/m000;
    	Iyz=m011-m010*m001/m000;

	// area-weighted moment of inertia
	cm << m100/m000,m010/m000,m001/m000;

    	double matrix[3][3] = {{Ixx,-Ixy,-Ixz},{-Ixy,Iyy,-Iyz},{-Ixz,-Iyz,Izz}};
    	double evecMat[3][3];
    	double eVal[3];
    	int minPos;

    	eigen_decomposition(matrix,evecMat,eVal);

	Vector3d radii;
	for(int i = 0; i< 3; i++)
	radii(i,0) = 1.0/sqrt(fabs(eVal[i]));

	// sort principle axis
	vector<float> orig;
	for(int i=0;i<3;i++)
	orig.push_back(radii(i,0));

	// this is a vector of {value,index} pairs
	vector<pair<float,size_t> > vp;
	vp.reserve(orig.size());

	for (size_t i = 0 ; i != orig.size() ; i++)
    	vp.push_back(make_pair(orig[i], i));

	// sorting will put lower values ahead of larger ones, resolving ties using the original index
	sort(vp.begin(), vp.end());

	for(int i=0 ; i<3 ; i++)
	{
    		pcx(i,0) = evecMat[i][vp[0].second];
		pbx(i,0) = evecMat[i][vp[1].second];
                pax(i,0) = evecMat[i][vp[2].second];
	}

	// change the global coordiantes w.r.t. centroid
    	for(int gno=0;gno<Gvertex.size();gno++)
    	{
        	Gvertex[gno]->x-=cm(0,0);
        	Gvertex[gno]->y-=cm(1,0);
        	Gvertex[gno]->z-=cm(2,0);
    	}

    	// make the surface normals (to each triangle elements) oriented outwards
	double dsum(0.0);
        for(int i =0;i< 3;i++)
	{
                x = Gvertex[element3D[0]->side[i].orientation[1]]->x - Gvertex[element3D[0]->side[i].orientation[0]]->x;
                y = Gvertex[element3D[0]->side[i].orientation[1]]->y + Gvertex[element3D[0]->side[i].orientation[0]]->y;
                dsum += x*y;
	}

	// sign of the area of the triangle is negative (<0), i.e. triangle vertex orientation clock-wise, normal needs to be flipped
        if(dsum > 0)
	{
                for(int eno=0;eno<element3D.size();eno++)
        	element3D[eno]->givenNormal*=-1.0;

	}
}

void edgeDescriptors(vector<Vertics*>& Gvertex,vector<elementList*> & eList,vector<Triangle3D*>& element3D)
{
        // edge descriptors

	int eno,neno,p1,p2;
    	double x,y,z,signal,edgeAngle,cosTheta,sinTheta;

	// calculate the 3D-edge angles and angle-axis pairs
    	for(int i=0;i<eList.size();i++)
    	{
        	x = Gvertex[eList[i]->nodes[0]]->x-Gvertex[eList[i]->nodes[1]]->x;
        	y = Gvertex[eList[i]->nodes[0]]->y-Gvertex[eList[i]->nodes[1]]->y;
        	z = Gvertex[eList[i]->nodes[0]]->z-Gvertex[eList[i]->nodes[1]]->z;

                // length of an edge
        	eList[i]->length = sqrt(x*x+y*y+z*z);

        	eno = eList[i]->sharedElement[0];
        	neno = eList[i]->sharedElement[1];


       	 	p2 = element3D[eno]->side[eList[i]->sharedEdge[0]].excludePoint;
        	p1 = element3D[neno]->side[eList[i]->sharedEdge[1]].excludePoint;

        	x = Gvertex[p1]->x-Gvertex[p2]->x;
        	y = Gvertex[p1]->y-Gvertex[p2]->y;
        	z = Gvertex[p1]->z-Gvertex[p2]->z;

        	Vector3d e(x,y,z);
        	e/=e.norm();

                // this demands the triangle normal MUST BE OUTWARD !!!
        	signal = element3D[eno]->givenNormal.dot(e);

		double acosCorr = element3D[eno]->givenNormal.dot(element3D[neno]->givenNormal);
		edgeAngle = acos(min(max(acosCorr,-1.0),1.0));

		// assign edge bending angle
                element3D[eno]->side[eList[i]->sharedEdge[0]].edgBendAngle = edgeAngle;
		element3D[neno]->side[eList[i]->sharedEdge[1]].edgBendAngle = edgeAngle;

        	// edge angle convex
        	if(signal < 0)
        	edgeAngle = PI - edgeAngle;

                // edge angle concave
        	else
        	edgeAngle = PI + edgeAngle;

		eList[i]->edge3DAngle = edgeAngle;

		// assign triangle edge-3d-angle (PI +/- edge-bending-angle)
                element3D[eno]->side[eList[i]->sharedEdge[0]].edgAngle = eList[i]->edge3DAngle;
                element3D[neno]->side[eList[i]->sharedEdge[1]].edgAngle = eList[i]->edge3DAngle;
    	}

	/// Calculate perimeter for each triangle elements
    	for(int eno=0;eno< element3D.size(); eno++)
    	{
        	double l(0.0);

        	for(int sid=0;sid<3;sid++)
        	{
            		l+= eList[element3D[eno]->sideToEdge[sid]]->length;

                        double sidePcat(element3D[eno]->side[sid].edgBendAngle/element3D[eno]->edgAngNorm);

			element3D[eno]->side[sid].pCat = sidePcat;
        	}
        	element3D[eno]->periMeter=l;
   	 }
}

void connectWithGlobe(vector<Vertics*>& Gvertex,vector<Triangle3D*>& element3D)
{
	// link the global coordinates to the triangles
	for(int eno=0;eno<element3D.size();eno++)
    	{
        	element3D[eno]->Vertex[0] << Gvertex[element3D[eno]->vertexIds[0]]->x,Gvertex[element3D[eno]->vertexIds[0]]->y,Gvertex[element3D[eno]->vertexIds[0]]->z;
        	element3D[eno]->Vertex[1] << Gvertex[element3D[eno]->vertexIds[1]]->x,Gvertex[element3D[eno]->vertexIds[1]]->y,Gvertex[element3D[eno]->vertexIds[1]]->z;
        	element3D[eno]->Vertex[2] << Gvertex[element3D[eno]->vertexIds[2]]->x,Gvertex[element3D[eno]->vertexIds[2]]->y,Gvertex[element3D[eno]->vertexIds[2]]->z;

                // mid-point coordinate for each triangle
		element3D[eno]->midPoint = (element3D[eno]->Vertex[0]+element3D[eno]->Vertex[1]+element3D[eno]->Vertex[2])/3.0;
    	}
}

void Image3dTo2D(vector<Region*> &regions,vector<Triangle3D*> &element3D, double pCatCurvList[])
{
    // rotate all the 3d-triangles and put them on x-y plane (normal z-axis)
    Vector3d z;
    z << 0.0,0.0,1.0;

    Vector3d cross;
    double cosAngle(0.0);
    double crossMagnitude(0.0);
    bool match(false);

    for(int eno = 0; eno < element3D.size(); eno++)
    {
        // quarternion rotation angle (cosine)
        cosAngle = element3D[eno]->givenNormal.dot(z);

        // quaternion rotation axis
        cross = element3D[eno]->givenNormal.cross(z);
        crossMagnitude = cross.norm();

        Vector3d nz(0.0,0.0,0.0);

        // assign quarternion rotation angle to this triangle region
        regions[eno]->rotAngle = acos(cosAngle);

        // the 3d-triangle required a rotation (triangle normal is non-parallel to z-axis)
        if((crossMagnitude > 0) || (cosAngle < 0))
        {
            // triangle normal is not anti-parallel to z-axis
            if(crossMagnitude > 0)
            {
                // nomalize quarternion rotation axis
                cross/=crossMagnitude;
                cross*=sin(regions[eno]->rotAngle/2);
            }

            // triangle normal is anti-parallel to z-axis
            else
            cross<< 0.0,-1.0,0.0;

            // get the quarternion of this triangle region
            regions[eno]->Q.w()=cos(regions[eno]->rotAngle/2);
            regions[eno]->Q.x()=cross(0,0);
            regions[eno]->Q.y()=cross(1,0);
            regions[eno]->Q.z()=cross(2,0);

            // find the necessary Quaternion for 3d -> 2d triangle rotation and rotate the 3d-triangle to 2d-triangle
            Quaternion<double> qr(0.0,0.0,0.0,0.0);

            // rotate the mid-point of the triangle region
            Quaternion<double> qMp(0,element3D[eno]->midPoint(0,0),element3D[eno]->midPoint(1,0),element3D[eno]->midPoint(2,0));
            qr = regions[eno]->Q * qMp * regions[eno]->Q.inverse();
            regions[eno]->zOffset = qr.z();

            // rotate all other three (=3) vertices of the triangle region
            Quaternion<double> qA(0,element3D[eno]->Vertex[0](0,0),element3D[eno]->Vertex[0](1,0),element3D[eno]->Vertex[0](2,0));
            qr = regions[eno]->Q * qA * regions[eno]->Q.inverse();
            regions[eno]->Vertics[0](0,0)=qr.x();
            regions[eno]->Vertics[0](1,0)=qr.y();

            Quaternion<double> qB(0,element3D[eno]->Vertex[1](0,0),element3D[eno]->Vertex[1](1,0),element3D[eno]->Vertex[1](2,0));
            qr = regions[eno]->Q * qB * regions[eno]->Q.inverse();
            regions[eno]->Vertics[1](0,0)=qr.x();
            regions[eno]->Vertics[1](1,0)=qr.y();

            Quaternion<double> qC(0,element3D[eno]->Vertex[2](0,0),element3D[eno]->Vertex[2](1,0),element3D[eno]->Vertex[2](2,0));
            qr = regions[eno]->Q * qC * regions[eno]->Q.inverse();
            regions[eno]->Vertics[2](0,0)=qr.x();
            regions[eno]->Vertics[2](1,0)=qr.y();

	    // rotate the triangle normal it self (to check the rotation accuracy: should be parallel to z-axis )
            Quaternion<double> nC(0,element3D[eno]->givenNormal(0,0),element3D[eno]->givenNormal(1,0),element3D[eno]->givenNormal(2,0));
            qr = regions[eno]->Q * nC * regions[eno]->Q.inverse();
            nz << qr.x(),qr.y(),qr.z();
        }

        // there is no need to rotate the 3d-triangle (triangle normal is parallel to z-axis)
        else
        {
            regions[eno]->Q.w()=0.0;
            regions[eno]->Q.x()=0.0;
            regions[eno]->Q.y()=0.0;
            regions[eno]->Q.z()=0.0;

            regions[eno]->Vertics[0](0,0)=element3D[eno]->Vertex[0](0,0);
            regions[eno]->Vertics[0](1,0)=element3D[eno]->Vertex[0](1,0);
            regions[eno]->Vertics[1](0,0)=element3D[eno]->Vertex[1](0,0);
            regions[eno]->Vertics[1](1,0)=element3D[eno]->Vertex[1](1,0);
            regions[eno]->Vertics[2](0,0)=element3D[eno]->Vertex[2](0,0);
            regions[eno]->Vertics[2](1,0)=element3D[eno]->Vertex[2](1,0);

            regions[eno]->zOffset = element3D[eno]->midPoint(2,0);
            nz = element3D[eno]->givenNormal;
        }

        // check the accuracy of the Quaternion rotation
        nz/=nz.norm();
        double s = nz.dot(z);

        if(s>0)
        match = true;

        else
        {
            match = false;
            cout << "Normals not aligned towards z (anti-parallel to) axis !!! emergency exit"<<endl;
            exit(-1);
        }

        // change the local coordinate frame to the Center of Mass frame
        regions[eno]->midPoint = (regions[eno]->Vertics[0]+regions[eno]->Vertics[1]+regions[eno]->Vertics[2])/3.0;
        for(int i=0;i<3;i++)
        regions[eno]->Vertics[i]-=regions[eno]->midPoint;

        // assign the necessary components from 3d-triangle to 2d-triangle
        regions[eno]->regionId = element3D[eno]->elementId;
        regions[eno]->side = element3D[eno]->side;
        regions[eno]->sideMap = element3D[eno]->sideMap;
        regions[eno]->sideRevMap = element3D[eno]->sideRevMap;
        regions[eno]->area = element3D[eno]->area;
        regions[eno]->periMeter = element3D[eno]->periMeter;
	regions[eno]->faceTag = element3D[eno]->faceTag;

        // assign appropriate pCat at the edges of the triangular region
	double totalPcat(element3D[eno]->pcatEdgeWeight*pCatCurvList[element3D[eno]->edgeTag]);


        // keep all the probabilities conserved (i.e in the range: 0.0 - 1.0)
        if(totalPcat > 1.0)
        totalPcat = 1.0;

	for(int sid=0;sid<3;sid++)
	regions[eno]->side[sid].pCat = totalPcat*element3D[eno]->side[sid].pCat;
    }
    Matrix2d A1,A2;
    Vector2d b1,b2;
    int v0,v1,v2;
    double cosTheta,sinTheta,x,y;
    match = false;

    // calculate the rotation angle between two 2d-triangle side
    for(int eno=0;eno<element3D.size();eno++)
    {
        for(int j=0;j<3;j++)
        {
            int ceno = regions[eno]->sideMap[j];
            int cesd = regions[ceno]->sideRevMap[regions[eno]->regionId];

            v1 = regions[eno]->side[j].dir[0];
            v2 = regions[eno]->side[j].dir[1];
            regions[eno]->side[j].midPoint  = (regions[eno]->Vertics[v1]+regions[eno]->Vertics[v2])/2.0;

            v0 = regions[eno]->side[j].excludePointMap[regions[eno]->side[j].excludePoint];

            A1 << (regions[eno]->Vertics[v1](0,0)-regions[eno]->Vertics[v0](0,0)),(regions[eno]->Vertics[v2](0,0)-regions[eno]->Vertics[v0](0,0)),
                  (regions[eno]->Vertics[v1](1,0)-regions[eno]->Vertics[v0](1,0)),(regions[eno]->Vertics[v2](1,0)-regions[eno]->Vertics[v0](1,0));

            b1 << regions[eno]->Vertics[v0](0,0),regions[eno]->Vertics[v0](1,0);

            Vector2d AB = regions[eno]->Vertics[v1]-regions[eno]->Vertics[v2];
            AB/=AB.norm();

            regions[eno]->side[j].dir2D = AB;
            v1 = regions[ceno]->side[cesd].dir[1];
            v2 = regions[ceno]->side[cesd].dir[0];

            v0 = regions[ceno]->side[cesd].excludePointMap[regions[ceno]->side[cesd].excludePoint];

            A2 << (regions[ceno]->Vertics[v1](0,0)-regions[ceno]->Vertics[v0](0,0)),(regions[ceno]->Vertics[v2](0,0)-regions[ceno]->Vertics[v0](0,0)),
                  (regions[ceno]->Vertics[v1](1,0)-regions[ceno]->Vertics[v0](1,0)),(regions[ceno]->Vertics[v2](1,0)-regions[ceno]->Vertics[v0](1,0));

            b2 << regions[ceno]->Vertics[v0](0,0),regions[ceno]->Vertics[v0](1,0);

            regions[eno]->side[j].A = A2*A1.inverse();
            regions[eno]->side[j].b = b2 - regions[eno]->side[j].A*b1;

            regions[ceno]->side[cesd].A = A1*A2.inverse();
            regions[ceno]->side[cesd].b = b1 - regions[ceno]->side[cesd].A*b2;

            MatrixXd G =A1.inverse();

            if(fabs(A1.determinant()*A2.determinant())>0.0)
            match = true;

            else
            {
                match = false;
                cout << "Affine matrix with Determinant = 0.0, is detected !!! emergency exit"<<endl;
                exit(-1);
            }


            Vector2d BA = regions[ceno]->Vertics[v1]-regions[ceno]->Vertics[v2];
            BA/=BA.norm();

            Vector3d s1(AB(0,0),AB(1,0),0.0),s2(BA(0,0),BA(1,0),0.0),ts(0.0,0.0,0.0);

            ts = s1.cross(s2);

            cosTheta = s1.dot(s2);
            sinTheta = ts.norm();

            Vector3d z(0.0,0.0,1.0);
            if(ts.norm()>0)
            {
                int k;
                ts/=ts.norm();
                k = (int) round(z.dot(ts));
                sinTheta*=k;
            }

            double angle = atan2(sinTheta,cosTheta);

            if(angle<0)
            angle += 2*PI;
            if(angle>2*PI)
            angle-=2*PI;

            regions[eno]->side[j].xyRotationAngle = angle;
        }

        // get the orientation axis (required for calculating the order parameter)
	vector<Vector3d> ends3D;
    	for(int p = 0;p < 2;p++)
    	ends3D.push_back(Vector3d(0,0,0));

        for(int i=0;i<3;i++)
        regions[eno]->orientation.push_back(Vector3d(0.0,0.0,0.0));

        x = 1.0;
        y = 0.0;

        // first two (=2) base vactor : rotate 2-d (x-axis and y-axis) -> triangle normal plane
        for(int ax = 0;ax < 2;ax++)
        {
            cross << regions[eno]->Q.x(),regions[eno]->Q.y(),regions[eno]->Q.z();
            crossMagnitude = cross.norm();

            // if rotation required
            if(crossMagnitude>0.0)
            {
                Quaternion<double> qr(0.0,0.0,0.0,0.0);

                // rotate the base vector from 2d-plane to 3d-triangle plane
                Quaternion<double> qI0(0,0.0,0.0,0.0);
                qr = regions[eno]->Q.inverse() * qI0 *regions[eno]->Q;
                ends3D[0](0,0)=qr.x();
                ends3D[0](1,0)=qr.y();
                ends3D[0](2,0)=qr.z();

                Quaternion<double> qI1(0,x,y,0.0);
                qr = regions[eno]->Q.inverse() * qI1 * regions[eno]->Q;
                ends3D[1](0,0)=qr.x();
                ends3D[1](1,0)=qr.y();
                ends3D[1](2,0)=qr.z();
            }

	    // if no rotation required
            else
            {
                ends3D[0](0,0)= 0.0;
                ends3D[0](1,0)= 0.0;
                ends3D[0](2,0)= 0.0;

                ends3D[1](0,0)= x;
                ends3D[1](1,0)= y;
                ends3D[1](2,0)= 0.0;
            }

            regions[eno]->orientation[ax]<< (ends3D[1](0,0)-ends3D[0](0,0)),(ends3D[1](1,0)-ends3D[0](1,0)),(ends3D[1](2,0)-ends3D[0](2,0));

            // change base vector from x-axis[(0,0),(1,0)] to y-axis[(0,0),(0,1)]
            swap(x,y);
        }

        // 3rd base vector (z-axis -> triangle normal map)
        regions[eno]->orientation[2]=element3D[eno]->givenNormal;

	ends3D.clear();
    }

}

void viewGraph(vector<Vertics*>& Gvertex,vector<Triangle3D*>& element3D,Vector3d &centroid,string outputDir)
{
    	string filename1(outputDir+"/ViewMeshDefault.off");
	string filename2(outputDir+"/ViewMeshFaceTag.dat");
	string filename3(outputDir+"/ViewMeshEdgeTag.dat");
	fstream fp1,fp2,fp3;

    	fp1.open(filename1.c_str(),ios::out);
        fp2.open(filename2.c_str(),ios::out);
	fp3.open(filename3.c_str(),ios::out);

	fp1 << "OFF" <<endl;
	fp1 << Gvertex.size() << " "<< element3D.size() << " " << 0 <<endl;

	for(int i=0;i<Gvertex.size();i++)
	fp1 << Gvertex[i]->x + centroid(0,0) << " " << Gvertex[i]->y + centroid(1,0) << " " << Gvertex[i]->z + centroid(2,0) <<endl;

        for(int i=0;i<element3D.size();i++)
	{
		fp1 << 3 << " ";
		for(int j =0;j<3;j++)
		fp1 << element3D[i]->vertexIds[j] << " ";
		fp1  <<endl;

		fp2 << element3D[i]->faceTag <<" ";
		fp3 << element3D[i]->edgeTag <<" ";
	}

    	fp1.close();
	fp2.close();
	fp3.close();
}

bool linePlaneIntersect(Vector3d n, Vector3d p, Vector3d a, Vector3d b, Vector3d& I)
{
    Vector3d ba = b-a;
    float nDotA = n.dot(a);
    float nDotBA = n.dot(ba);
    float d = n.dot(p);

    float t = (d - nDotA)/nDotBA;
    bool intersect(false);

    if (t >= 0.0f && t <= 1.0f)
    {
        intersect = true;
        I = a + t*ba;
    }

    return(intersect);
}

double areaPolygon3D(vector<Vertics>& V,Vector3d N)
{
    int n = V.size();
    V.push_back(Vertics(V[0].x,V[0].y,V[0].z));

    float area = 0;
    float an, ax, ay, az;
    int  coord;
    int  i, j, k;

    if (n < 3) return 0;

    // select largest abs coordinate to ignore for projection
    ax = (N(0,0)>0 ? N(0,0) : -N(0,0));
    ay = (N(1,0)>0 ? N(1,0) : -N(1,0));
    az = (N(2,0)>0 ? N(2,0) : -N(2,0));

    coord = 3;
    if (ax > ay)
    {
        if (ax > az)
        coord = 1;
    }
    else if (ay > az) coord = 2;

    // compute area of the 2D projection
    switch (coord)
    {
      case 1:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].y * (V[j].z - V[k].z));
        break;

      case 2:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].z * (V[j].x - V[k].x));
        break;

      case 3:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].x * (V[j].y - V[k].y));
        break;
    }

    switch (coord)
    {
      case 1:
        area += (V[n].y * (V[1].z - V[n-1].z));
        break;
      case 2:
        area += (V[n].z * (V[1].x - V[n-1].x));
        break;
      case 3:
        area += (V[n].x * (V[1].y - V[n-1].y));
        break;
    }

    // scale to get area before projection
    an = sqrt( ax*ax + ay*ay + az*az);
    switch (coord)
    {
      case 1:
        area *= (an / (2 * N(0,0)));
        break;
      case 2:
        area *= (an / (2 * N(1,0)));
        break;
      case 3:
        area *= (an / (2 * N(2,0)));
    }

    return fabs(area);
}

double intersectingPolygon(vector<Vertics>& polygon,vector<Vertics>& Gvertex,vector<Region*>& regions,vector<elementList> & eList,Vector3d np, Vector3d cp, string shapeType)
{
	double polyArea(0.0);
        if(shapeType=="2D-plane")
	{
		cout << "2d-plane, so no division plane" <<endl;
	}
	else
	{
		polygon.clear();

		for(int eno=0;eno< regions.size();eno++)
		regions[eno]->polyIntersectMark==0;

		vector<int> intersectedEdges;
		double x,y,z;
		for(int i =0;i< eList.size();i++)
		{
		        x = Gvertex[eList[i].nodes[0]].x;
		        y = Gvertex[eList[i].nodes[0]].y;
			z = Gvertex[eList[i].nodes[0]].z;
		        Vector3d lS(x,y,z);

			x = Gvertex[eList[i].nodes[1]].x;
		        y = Gvertex[eList[i].nodes[1]].y;
			z = Gvertex[eList[i].nodes[1]].z;
		        Vector3d lE(x,y,z);

			Vector3d I(0.0,0.0,0.0);

			bool Intersect = linePlaneIntersect(np,cp,lS,lE,I);

			if(Intersect)
		        {
		       		intersectedEdges.push_back(i);
				eList[i].intersectionByPlane = I;
		        }
		}

		for(int i =0;i< intersectedEdges.size();i++)
		{
		        int edg = intersectedEdges[i];

		        int faceID = eList[edg].sharedElement[0];
		        regions[faceID]->intersectEdg.push_back(edg);
		        regions[faceID]->polyIntersectMark = 1;

		        faceID = eList[edg].sharedElement[1];
		        regions[faceID]->intersectEdg.push_back(edg);
		        regions[faceID]->polyIntersectMark = 1;
		}

		vector<int> polygonEdg;
		int startEdgeIndex = intersectedEdges[0];
		int currentEdgeIndex = startEdgeIndex;

		polygonEdg.push_back(currentEdgeIndex);
		int currentFace = eList[currentEdgeIndex].sharedElement[0];

		while (true)
		{
			currentEdgeIndex = (currentEdgeIndex == regions[currentFace]->intersectEdg[0]) ?  regions[currentFace]->intersectEdg[1] : regions[currentFace]->intersectEdg[0];

			currentFace = (currentFace == eList[currentEdgeIndex].sharedElement[0]) ?  eList[currentEdgeIndex].sharedElement[1] : eList[currentEdgeIndex].sharedElement[0];

			if(currentEdgeIndex == startEdgeIndex)
		    	break;

		        polygonEdg.push_back(currentEdgeIndex);
		}

		// construct polygon
		for(int i =0; i< polygonEdg.size();i++)
		{
		        int edg = polygonEdg[i];
			Vector3d pv = eList[edg].intersectionByPlane;
			polygon.push_back(Vertics(pv(0,0),pv(1,0),pv(2,0)));
		}

		// polygon area
		polyArea = areaPolygon3D(polygon,np);

		intersectedEdges.clear();
		polygonEdg.clear();

		for(int i=0;i<  regions.size() ;i++)
		regions[i]->intersectEdg.clear();

		return(polyArea);
	}
}

vector<string> split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;

  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }

  return internal;
}

/// PBC for 2D-plane and disc
void establishPBC(vector<Vertics*>& vertices,vector<Triangle3D*>& triangles,vector<Triangle3D*>& viewTriangles,vector<elementList*>& eList,double sizeFactor)
{
	int bCount(0);
	vector<int>edgeInMesh,edgeInBoundary;

	/// find the boundary edges
	for(int i=0;i<eList.size();i++)
	{
		bool hitorigin(eList[i]->nodes[1]==5);

		if(hitorigin)
		bCount++;

		if(!hitorigin)
		edgeInMesh.push_back(i);
	}

	int firstBT=triangles.size()-bCount;

	for(int i=0;i<edgeInMesh.size();i++)
	{
		int n1 = eList[edgeInMesh[i]]->sharedElement[0];
		int n2 = eList[edgeInMesh[i]]->sharedElement[1];
		int largest = ((n1 > n2) ? n1 : n2);

		bool inBoundary(largest>=firstBT);

		if(inBoundary)
		edgeInBoundary.push_back(edgeInMesh[i]);
	}

	edgeInMesh.clear();

	/// delete PPB communicating triangles
	for(int i=0;i<bCount;i++)
	{
		triangles.pop_back();
		viewTriangles.pop_back();
	}

	/// boundary edge (clock/anti-clock wise) sorting
	int vid(vertices.size()-1);
	for(int i=0;i<edgeInBoundary.size();i++)
	{
		vid++;
		eList[edgeInBoundary[i]]->midPoint = vid;

		double xm = (vertices[eList[edgeInBoundary[i]]->nodes[0]]->x+vertices[eList[edgeInBoundary[i]]->nodes[1]]->x)/2.0;
		double ym = (vertices[eList[edgeInBoundary[i]]->nodes[0]]->y+vertices[eList[edgeInBoundary[i]]->nodes[1]]->y)/2.0;
		double zm = (vertices[eList[edgeInBoundary[i]]->nodes[0]]->z+vertices[eList[edgeInBoundary[i]]->nodes[1]]->z)/2.0;

		Vector3d vm(xm,ym,zm);
		vertices.push_back(new Vertics(vm(0,0),vm(1,0),vm(2,0)));
	}

	double elementDistance;
	vector<double> edgePositionBoundary;

	for(int i=0;i<edgeInBoundary.size();i++)
	{
                Vector2d radialExt(vertices[eList[edgeInBoundary[i]]->midPoint]->x,vertices[eList[edgeInBoundary[i]]->midPoint]->y);
                radialExt/=radialExt.norm();

		double theta = atan2(radialExt(1,0),radialExt(0,0));

		if(theta<0)
		theta+=2*PI;

		if(theta>2*PI)
		theta-=2*PI;

		edgePositionBoundary.push_back(theta);
	}

	/// remove un-necessary vertices
        for(int i=0;i<=edgeInBoundary.size();i++)
        vertices.pop_back();

	vector<pair<double,int> > vp;
	vp.reserve(edgePositionBoundary.size());

	for(int i = 0 ; i != edgePositionBoundary.size() ; i++)
	vp.push_back(make_pair(edgePositionBoundary[i],i));
	edgePositionBoundary.clear();

	sort(vp.begin(),vp.end());

	vector<int> indexBoundary;
	for(int i=0;i<vp.size();i++)
	indexBoundary.push_back(vp[i].second);

	int n1,n2,p1,p2,e1,e2,s1,s2,l1,l2;

	/// establish top-botton PBC connectivity
	e1= edgeInBoundary[indexBoundary[0]];
	n1 = eList[e1]->sharedElement[0];
	n2 = eList[e1]->sharedElement[1];
	s1 = eList[e1]->sharedEdge[0];
	s2 = eList[e1]->sharedEdge[1];
	if(n1<n2)
	{
		p1=n1;
		l1 = s1;
	}
	else
	{
		p1=n2;
		l1 = s2;
	}
	e2 = edgeInBoundary[indexBoundary[2]];
	n1 = eList[e2]->sharedElement[0];
	n2 = eList[e2]->sharedElement[1];
	s1 = eList[e2]->sharedEdge[0];
	s2 = eList[e2]->sharedEdge[1];
	if(n1<n2)
	{
		p2=n1;
		l2 = s1;
	}
	else
	{
		p2=n2;
		l2 = s2;
	}

	triangles[p1]->periodicTag = p2;
	triangles[p2]->periodicTag = p1;
	triangles[p1]->sideMap[l1] = p2;
	triangles[p2]->sideMap[l2] = p1;
	triangles[p1]->sideRevMap[p2] = l1;
	triangles[p2]->sideRevMap[p1] = l2;
	triangles[p1]->sideToEdge[l1] = e1;
	triangles[p2]->sideToEdge[l2] = e2;

	/// establish left-right PBC connectivity
	e1= edgeInBoundary[indexBoundary[1]];
	n1 = eList[e1]->sharedElement[0];
	n2 = eList[e1]->sharedElement[1];
	s1 = eList[e1]->sharedEdge[0];
	s2 = eList[e1]->sharedEdge[1];
	if(n1<n2)
	{
		p1=n1;
		l1 = s1;
	}
	else
	{
		p1=n2;
		l1 = s2;
	}
	e2 = edgeInBoundary[indexBoundary[3]];
	n1 = eList[e2]->sharedElement[0];
	n2 = eList[e2]->sharedElement[1];
	s1 = eList[e2]->sharedEdge[0];
	s2 = eList[e2]->sharedEdge[1];
	if(n1<n2)
	{
		p2=n1;
		l2 = s1;
	}
	else
	{
		p2=n2;
		l2 = s2;
	}
	triangles[p1]->periodicTag = p2;
	triangles[p2]->periodicTag = p1;
	triangles[p1]->sideMap[l1] = p2;
	triangles[p2]->sideMap[l2] = p1;
	triangles[p1]->sideRevMap[p2] = l1;
	triangles[p2]->sideRevMap[p1] = l2;
	triangles[p1]->sideToEdge[l1] = e1;
	triangles[p2]->sideToEdge[l2] = e2;

    	indexBoundary.clear();
}

