#include "mesh.hpp"
#include "geometry.hpp"
#include "region.hpp"
#include "system.hpp"
#include "mt_tip.hpp"
#include "linalg.hpp"

// initialize edge
void iniEdgeRecord(struct EdgeRecord* head, vector<int>& v, int cid, int fid)
{

    for (int i = 0; i < 2; i++)
    {
        head->key.push_back(v[i]);
    }

    head->cell_id = cid;
    head->edge_index = fid;
    head->next = nullptr;
}

// appending edge
void addEdgeRecord(struct EdgeRecord* head, vector<int>& v, int cid, int fid)
{
    EdgeRecord* newedgeRecord = new EdgeRecord;
    newedgeRecord->key = v;
    newedgeRecord->cell_id = cid;
    newedgeRecord->edge_index = fid;
    newedgeRecord->next = nullptr;
    EdgeRecord* cur = head;

    while (cur)
    {
        if (cur->next == nullptr)
        {
            cur->next = newedgeRecord;
            return;
        }
        cur = cur->next;
    }
}

// delete from edge front
void deleteFirst(struct EdgeRecord** head)
{
    struct EdgeRecord* tmp = *head;
    if (tmp == nullptr)
    {
        return;
    }
    *head = tmp->next;
    delete tmp;
    tmp = nullptr;
}

// display edge record
void display(struct EdgeRecord* head)
{
    EdgeRecord* list = head;
    int listCount(0);
    while (list)
    {
        cout << "-----------------------------" << endl;
        cout << "Member: " << listCount << endl;
        cout << "-----------------------------" << endl;

        cout << "The key ids: ";
        for (int i = 0; i < 2; i++)
        {
            cout << list->key[i] << " ";
        }

        cout << endl;
        cout << "The cell_id: " << list->cell_id << endl;
        cout << "The face_index: " << list->edge_index << endl;
        list = list->next;
        listCount++;
    }
    cout << endl;
}

// bin-sorting
void edge_sort(EdgeRecord** e_list, int node_max)
{
    vector<int> v1;
    for (int i = 0; i < 2; i++)
    {
        v1.push_back(0);
    }

    for (int i = 1; i >= 0; i--)
    {
        EdgeRecord** Bin;
        Bin = new EdgeRecord*[node_max];

        for (int j = 0; j < node_max; j++)
        {
            Bin[j] = new EdgeRecord;
            iniEdgeRecord(Bin[j], v1, -1, -1);
        }

        EdgeRecord* edge = *e_list;
        while (edge)
        {
            addEdgeRecord(Bin[edge->key[i]], edge->key, edge->cell_id, edge->edge_index);

            EdgeRecord* oldedge = edge;
            edge = edge->next;
            delete oldedge;
            oldedge = nullptr;
        }

        for (int j = 0; j < node_max; j++)
        {
            deleteFirst(&Bin[j]);
        }

        *e_list = new EdgeRecord;
        iniEdgeRecord(*e_list, v1, -1, -1);

        for (int j = 0; j < node_max; j++)
        {
            EdgeRecord* binedge = Bin[j];
            while (binedge)
            {
                addEdgeRecord(*e_list, binedge->key, binedge->cell_id, binedge->edge_index);
                binedge = binedge->next;
            }
        }
        deleteFirst(&(*e_list));
        delete[] Bin;
        Bin = nullptr;
    }
}

// connectivity establishment
void find_connectivity(EdgeRecord** f_list, vector<ElementList*>& e_list)
{
    EdgeRecord* edge = *f_list;
    int ecount(0);
    while (edge != nullptr)
    {
        EdgeRecord* edgenext = edge->next;
        if (edgenext == nullptr)
        {
            break;
        }

        if (edge->key == edgenext->key)
        {
            e_list.push_back(new ElementList());
            ecount = e_list.size() - 1;
            e_list[ecount]->nodes[0] = edge->key[0];
            e_list[ecount]->nodes[1] = edge->key[1];
            e_list[ecount]->sharedElement[0] = edge->cell_id;
            e_list[ecount]->sharedElement[1] = edgenext->cell_id;
            e_list[ecount]->sharedEdge[0] = edge->edge_index;
            e_list[ecount]->sharedEdge[1] = edgenext->edge_index;

            edge = edgenext->next;
        }
        else
        {
            edge = edge->next;
        }
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
    vector<Vertices*> Gvertex;

    // all triangle elemet list
    vector<Triangle3D*> element3D;

    // centroid/principle-axis of the cell
    Eigen::Vector3d centroid(0.0, 0.0, 0.0);
    Eigen::Vector3d principleAxis(0.0, 0.0, 0.0);
    Eigen::Vector3d axB(0.0, 0.0, 0.0);
    Eigen::Vector3d axC(0.0, 0.0, 0.0);

    // list of triangle-edges to be detected
    vector<ElementList*> eList;

    int a0(0), a1(0), a2(0);
    int vmax(0), trimax(0), spf(0), cur(0);
    string str;

    string filename;
    double x(0), y(0), z(0), eNorm, eAng(0.0);

    if (inFile == "cubeReal")
    {
        Gvertex.push_back(new Vertices(0.5 * sizeFactor, 0.5 * sizeFactor, 1.0 * sizeFactor)); // 0
        Gvertex.push_back(new Vertices(0.5 * sizeFactor, 0.5 * sizeFactor, 0.0));              // 1

        Gvertex.push_back(new Vertices(0.5 * sizeFactor, 0.0, 0.5 * sizeFactor));              // 2
        Gvertex.push_back(new Vertices(0.5 * sizeFactor, 1.0 * sizeFactor, 0.5 * sizeFactor)); // 3

        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 0.5 * sizeFactor, 0.5 * sizeFactor)); // 4
        Gvertex.push_back(new Vertices(0.0, 0.5 * sizeFactor, 0.5 * sizeFactor));              // 5

        Gvertex.push_back(new Vertices(0.0, 0.0, 1.0 * sizeFactor));                           // 6
        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 0.0, 1.0 * sizeFactor));              // 7
        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 1.0 * sizeFactor, 1.0 * sizeFactor)); // 8
        Gvertex.push_back(new Vertices(0.0, 1.0 * sizeFactor, 1.0 * sizeFactor));              // 9

        Gvertex.push_back(new Vertices(0.0, 0.0, 0.0));                           // 10
        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 0.0, 0.0));              // 11
        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 1.0 * sizeFactor, 0.0)); // 12
        Gvertex.push_back(new Vertices(0.0, 1.0 * sizeFactor, 0.0));              // 13

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

        for (size_t i = 0; i < element3D.size(); i++)
        {
            element3D[i]->elementId = i;
        }
    }

    else if (inFile == "2D-plane")
    {
        Gvertex.push_back(new Vertices(0.0, 0.0, 0.0));
        Gvertex.push_back(new Vertices(0.0, 1.0 * sizeFactor, 0.0));
        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 1.0 * sizeFactor, 0.0));
        Gvertex.push_back(new Vertices(1.0 * sizeFactor, 0.0, 0.0));
        Gvertex.push_back(new Vertices(0.5 * sizeFactor, 0.5 * sizeFactor, 0.0));
        Gvertex.push_back(new Vertices(0.5 * sizeFactor, 0.5 * sizeFactor, -1.0));

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
        string inputMeshFile = g->system->p.inputDir + "/" + g->system->p.geometry + "/" + inFile + ".dat";

        // open input mesh file
        fstream fp;
        fp.open(inputMeshFile.c_str(), ios::in);

        fp >> vmax >> trimax;

        // read the global vertex coordinates
        for (int i = 0; i < vmax; i++)
        {
            fp >> x >> y >> z;
            Gvertex.push_back(new Vertices(0.0, 0.0, 0.0));

            // scale the mesh, if asked
            Gvertex[i]->x = x * sizeFactor;
            Gvertex[i]->y = y * sizeFactor;
            Gvertex[i]->z = z * sizeFactor;
        }

        // read triangle informatio:
        // (a) vertex coordinates
        // (b) spcial face tag
        // (c) curvature tag
        // (d) Pcat norm factor
        // (d) average edge angle (of the curvature domain it belongs)
        for (int i = 0; i < trimax; i++)
        {
            fp >> a0 >> a1 >> a2 >> spf >> cur >> eNorm >> eAng;
            element3D.push_back(new Triangle3D());

            element3D[i]->elementId = i;
            element3D[i]->faceTag = spf;
            element3D[i]->edgeTag = cur;
            element3D[i]->edgAngNorm = eNorm;
            element3D[i]->pcatEdgeWeight = eAng;

            // vertices need zero-based numbering
            element3D[i]->vertexIds[0] = a0 - 1;
            element3D[i]->vertexIds[1] = a1 - 1;
            element3D[i]->vertexIds[2] = a2 - 1;
        }

        fp >> str;

        // close the input mesh file
        fp.close();
    }

    // establish connectivity among the triangles
    makeGraph(Gvertex.size(), eList, element3D);

    // store the numer of triangles, vertices and edges
    int F = element3D.size();
    int V = Gvertex.size();
    int E = eList.size();

    // check mesh consistency
    bool convex(checkEular(F, V, E));

    // copy triangle elements for viewing mesh (visualization in meshLab)
    vector<Triangle3D*> viewTriangle;
    for (size_t i = 0; i < element3D.size(); i++)
    {
        viewTriangle.push_back(new Triangle3D());
        viewTriangle[i]->vertexIds[0] = element3D[i]->vertexIds[0];
        viewTriangle[i]->vertexIds[1] = element3D[i]->vertexIds[1];
        viewTriangle[i]->vertexIds[2] = element3D[i]->vertexIds[2];

        viewTriangle[i]->faceTag = element3D[i]->faceTag;
        viewTriangle[i]->edgeTag = element3D[i]->edgeTag;
    }

    // check surface orientation
    orientSurface(Gvertex, element3D, viewTriangle);

    // check surface orientation
    bool surfaceOriented(checkSurfaceOrientation(element3D));

    // measure rigidbody properties
    rigidBodyProperties(Gvertex, element3D, centroid, principleAxis, axB, axC);

    // assign coordinate values to the triangle vertices
    connectWithGlobe(Gvertex, element3D);

    // get triangle edge related measures/properties
    edgeDescriptors(Gvertex, eList, element3D);

    /// PBC for 2D-plane
    if (inFile == "2D-plane")
    {
        principleAxis << 1.0, 0.0, 0.0;
        axB << 0.0, 1.0, 0.0;
        axC << 0.0, 0.0, 0.0;
        establishPBC(Gvertex, element3D, viewTriangle, eList, sizeFactor);
    }

    // visualize mesh in meshLab
    viewGraph(Gvertex, viewTriangle, centroid, g->system->p.outputDir + "/" + g->system->p.geometry + "/outputData");

    double totalArea(0.0);
    for (size_t eno = 0; eno < element3D.size(); eno++)
    {
        totalArea += element3D[eno]->area;
    }

    // show mesh in enabled
    if (showMesh == 1)
    {
        cout << "============================================" << endl;

        // show principle axis
        cout << "Principle Axis ---> (" << principleAxis(0, 0) << "," << principleAxis(1, 0) << ","
             << principleAxis(2, 0) << ")" << endl;

        // show the numer of triangles, vertices and edges
        cout << "F = " << F << " V = " << V << " E = " << E << endl;

        // mesh is non-convex
        if (!convex)
        {
            cout << "The surface mesh is non-Euler" << endl;
        }

        // mesh in convex
        else
        {
            cout << "The surface mesh is Euler" << endl;
        }

        if (surfaceOriented)
        {
            cout << "Surface oriented (outwards) successfully" << endl;
        }

        cout << "============================================" << endl;
    }

    // summary of the mesh properties
    fstream fp;
    filename = g->system->p.outputDir + "/" + g->system->p.geometry + "/outputData" + string("/MeshInfo.txt");
    fp.open(filename.c_str(), ios::out);

    fp << "CellIndex\t" << g->system->p.geomParam << endl;
    fp << "TriangleNumbers\t" << element3D.size() << endl;
    fp << "Centroid\t" << centroid(0, 0) << "," << centroid(1, 0) << "," << centroid(2, 0) << endl;
    fp << "PA\t" << principleAxis(0, 0) << "," << principleAxis(1, 0) << "," << principleAxis(2, 0) << endl;
    fp << "PB\t" << axB(0, 0) << "," << axB(1, 0) << "," << axB(2, 0) << endl;
    fp << "PC\t" << axC(0, 0) << "," << axC(1, 0) << "," << axC(2, 0) << endl;
    fp << "Area\t" << totalArea << endl;

    fp.close();

    // starting simulation on the surface-mesh
    map<string, double>::iterator it;

    // list of edge catastrophe values
    int edgNumber = g->system->p.edgCatMap.size();
    double pCatList[edgNumber];
    for (int i = 0; i <= edgNumber; i++)
    {
        pCatList[i] = 0.0;
    }

    for (it = g->system->p.edgCatMap.begin(); it != g->system->p.edgCatMap.end(); it++)
    {
        vector<string> currEdg = split(it->first, '_');
        int edgIdx{ std::stoi(currEdg[1]) };
        double edgCatVal = it->second;
        pCatList[edgIdx] = edgCatVal;
    }

    // list of face catastrophe values
    int faceNumber = g->system->p.faceCatMap.size();
    for (int i = 0; i <= faceNumber; i++)
    {
        g->system->p.RegionKcatMultiplier.push_back(1.0);
    }

    for (it = g->system->p.faceCatMap.begin(); it != g->system->p.faceCatMap.end(); it++)
    {
        vector<string> currFace = split(it->first, '_');
        int faceIdx{ std::stoi(currFace[1]) };
        double faceCatVal = it->second;
        g->system->p.RegionKcatMultiplier[faceIdx] = faceCatVal;
    }

    // assign the centroid and principle axis to the geometry
    g->objectCM = centroid;
    g->objectPA = principleAxis;
    g->nucleousPosition << 0.0, 0.0, 0.0;

    // total surface area of the geometry
    g->area = totalArea;

    // create only a single path area (at first place)
    g->patchArea[0] = g->area;

    // include all triangles to a single path (at first place)
    for (size_t eno = 0; eno < element3D.size(); eno++)
    {
        g->RegionsIndex[0].push_back(eno);
    }

    // total number of region in the geometry
    g->elementMax = element3D.size();

    // to create PPB
    // (a)  global vertex to be used to create PPB
    for (size_t i = 0; i < Gvertex.size(); i++)
    {
        g->globalVertex.push_back(Vertices(0.0, 0.0, 0.0));
        g->globalVertex[i].x = Gvertex[i]->x;
        g->globalVertex[i].y = Gvertex[i]->y;
        g->globalVertex[i].z = Gvertex[i]->z;
    }

    // (b) element list
    for (size_t i = 0; i < eList.size(); i++)
    {
        g->ppbEdgeList.push_back(ElementList());

        for (size_t j = 0; j < 2; j++)
        {
            g->ppbEdgeList[i].nodes[j] = eList[i]->nodes[j];
            g->ppbEdgeList[i].sharedElement[j] = eList[i]->sharedElement[j];
            g->ppbEdgeList[i].sharedEdge[j] = eList[i]->sharedEdge[j];
        }
    }

    // assign memory for all the 2d-triangular regions of the geometry
    for (size_t eno = 0; eno < element3D.size(); eno++)
    {
        g->regions.push_back(new Triangle(element3D[eno]->area, g));
    }

    // copy 3d-triangle properties information to 2d-triangular regions
    Image3dTo2D(g->regions, element3D, pCatList);

    // 3d-triangle properties are transfered to 2d-triangular regions, so release all the memory occupied by 3d
    // triangles
    for (size_t eno = 0; eno < element3D.size(); eno++)
    {
        delete element3D[eno];
        // FOLLOW-UP
        // This used to be a comparison (element3D[eno] == NULL).
        // Check if the fix introduces any regressions.
        element3D[eno] = nullptr;
    }
    element3D.clear();

    // global vertices no longer needed, so release the memories occupied by the global vertices
    for (size_t gno = 0; gno < Gvertex.size(); gno++)
    {
        delete Gvertex[gno];
        Gvertex[gno] = nullptr;
    }
    Gvertex.clear();

    // edge list is no longer needed, so release the memories occupied by the edges
    for (size_t el = 0; el < eList.size(); el++)
    {
        delete eList[el];
        eList[el] = nullptr;
    }
    eList.clear();

    return;
}

void makeGraph(int gvertex, vector<ElementList*>& e_list, vector<Triangle3D*>& element3_d)
{
    // make the connectivity graph
    void edge_sort(EdgeRecord**, int);
    void find_connectivity(EdgeRecord**, vector<ElementList*>&);

    struct EdgeRecord* EdgList = new EdgeRecord();

    vector<int> v;
    v.push_back(0);
    v.push_back(0);

    iniEdgeRecord(EdgList, v, -1, -1);

    v.clear();

    for (size_t i = 0; i < element3_d.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            element3_d[i]->vertexIdsMap[element3_d[i]->vertexIds[j]] = j;
        }

        for (int j = 0; j < 3; j++)
        {
            v.push_back(element3_d[i]->vertexIds[j % 3]);
            v.push_back(element3_d[i]->vertexIds[(j + 1) % 3]);
            element3_d[i]->side[j].orientation = v;

            element3_d[i]->side[j].orientationMap[v[0]] = j % 3;
            element3_d[i]->side[j].orientationMap[v[1]] = (j + 1) % 3;
            element3_d[i]->side[j].dir.push_back(j % 3);
            element3_d[i]->side[j].dir.push_back((j + 1) % 3);

            element3_d[i]->side[j].excludePoint = element3_d[i]->vertexIds[(j + 2) % 3];
            element3_d[i]->side[j].excludePointMap[element3_d[i]->vertexIds[(j + 2) % 3]] = (j + 2) % 3;

            sort(v.begin(), v.end());

            addEdgeRecord(EdgList, v, i, j);
            v.clear();
        }
        element3_d[i]->elementId = i;
    }
    deleteFirst(&(EdgList));

    edge_sort(&EdgList, gvertex);
    find_connectivity(&EdgList, e_list);

    for (size_t i = 0; i < e_list.size(); i++)
    {
        element3_d[e_list[i]->sharedElement[0]]->sideMap[e_list[i]->sharedEdge[0]] = e_list[i]->sharedElement[1];
        element3_d[e_list[i]->sharedElement[1]]->sideMap[e_list[i]->sharedEdge[1]] = e_list[i]->sharedElement[0];

        element3_d[e_list[i]->sharedElement[0]]->sideRevMap[e_list[i]->sharedElement[1]] = e_list[i]->sharedEdge[0];
        element3_d[e_list[i]->sharedElement[1]]->sideRevMap[e_list[i]->sharedElement[0]] = e_list[i]->sharedEdge[1];

        element3_d[e_list[i]->sharedElement[0]]->sideToEdge[e_list[i]->sharedEdge[0]] = i;
        element3_d[e_list[i]->sharedElement[1]]->sideToEdge[e_list[i]->sharedEdge[1]] = i;
    }

    delete EdgList;
    EdgList = nullptr;

    return;
}

void orientSurface(vector<Vertices*>& gvertex, vector<Triangle3D*>& element3_d, vector<Triangle3D*>& view_triangle)
{
    // orient surface normals
    void correctOrientation(vector<Vertices*>&, Triangle3D*&, Triangle3D*&, int, int, int&);

    vector<Triangle3D*> TriangleAlive;
    TriangleAlive.push_back(element3_d[0]);

    vector<Triangle3D*> TriangleDead;

    while (!TriangleAlive.empty())
    {
        Triangle3D* nextTriangle = TriangleAlive.back();
        TriangleAlive.pop_back();

        TriangleDead.push_back(nextTriangle);

        for (int i = 0; i < 3; i++)
        {
            Triangle3D* neighTriangle = element3_d[nextTriangle->sideMap[i]];

            bool isPresent = find(TriangleDead.begin(), TriangleDead.end(), neighTriangle) != TriangleDead.end();
            if (!isPresent)
            {
                int flip(1);
                int n1 = neighTriangle->sideRevMap[nextTriangle->elementId];
                correctOrientation(gvertex, nextTriangle, neighTriangle, i, n1, flip);
                TriangleAlive.push_back(neighTriangle);
                if (flip == -1)
                {
                    view_triangle[nextTriangle->sideMap[i]]->vertexIds[1]
                    = element3_d[nextTriangle->sideMap[i]]->vertexIds[2];
                    view_triangle[nextTriangle->sideMap[i]]->vertexIds[2]
                    = element3_d[nextTriangle->sideMap[i]]->vertexIds[1];
                }
            }
        }
    }
}

void correctOrientation(
vector<Vertices*>& gvertex, Triangle3D*& next_triangle, Triangle3D*& neigh_triangle, int n1, int n2, int& flip)
{
    if (next_triangle->side[n1].orientation == neigh_triangle->side[n2].orientation)
    {
        for (int j = 0; j < 3; j++)
        {
            reverse(neigh_triangle->side[j].orientation.begin(), neigh_triangle->side[j].orientation.end());
            reverse(neigh_triangle->side[j].dir.begin(), neigh_triangle->side[j].dir.end());
            flip *= -1;
        }
    }

    return;
}

bool checkSurfaceOrientation(vector<Triangle3D*>& element3_d)
{
    bool match(false);
    for (size_t i = 0; i < element3_d.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            match = false;
            int n = element3_d[i]->sideMap[j];
            int s = element3_d[n]->sideRevMap[i];
            vector<int> v;
            v.push_back(element3_d[n]->side[s].orientation[1]);
            v.push_back(element3_d[n]->side[s].orientation[0]);

            if (element3_d[i]->side[j].orientation == v)
            {
                match = true;
            }

            if (!match)
            {
                cout << "Surface not oriented ... emergency exit !!!" << endl;
                exit(1);
            }
        }
    }

    return (match);
}

bool checkEular(int f, int v, int e)
{
    bool graph(((3 * f / 2) - e) == 0);

    if (!graph)
    {
        cout << "Inconsistency detected between the nuber of TRIANGLES and EDGES  ... emergency exit !!!" << endl;
        exit(-1);
    }

    graph = ((f + v - e) == 2);
    if (!graph)
    {
        cout << "Euler condition violated ... emergency exit !!!" << endl;
        exit(-1);
    }

    return (graph);
}

void rigidBodyProperties(vector<Vertices*>& gvertex,
                         vector<Triangle3D*>& element3_d,
                         Eigen::Vector3d& cm,
                         Eigen::Vector3d& pax,
                         Eigen::Vector3d& pbx,
                         Eigen::Vector3d& pcx)
{
    // rigidBodyProperties
    vector<Eigen::Vector3d> xRow1, yRow1, zRow1, xRow2, yRow2, zRow2, xRow3, yRow3, zRow3;
    double x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double x_2, y_2, z_2, xy, xz, yz;
    double x_3, y_3, z_3;
    double x_2y, x_2z, y_2x, y_2z, z_2y, z_2x, xyz;

    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        // calculate normals to each triangle elements
        x = gvertex[element3_d[eno]->side[0].orientation[1]]->x - gvertex[element3_d[eno]->side[0].orientation[0]]->x;
        y = gvertex[element3_d[eno]->side[0].orientation[1]]->y - gvertex[element3_d[eno]->side[0].orientation[0]]->y;
        z = gvertex[element3_d[eno]->side[0].orientation[1]]->z - gvertex[element3_d[eno]->side[0].orientation[0]]->z;
        Eigen::Vector3d RL(x, y, z);

        x = gvertex[element3_d[eno]->side[0].orientation[0]]->x - gvertex[element3_d[eno]->side[0].excludePoint]->x;
        y = gvertex[element3_d[eno]->side[0].orientation[0]]->y - gvertex[element3_d[eno]->side[0].excludePoint]->y;
        z = gvertex[element3_d[eno]->side[0].orientation[0]]->z - gvertex[element3_d[eno]->side[0].excludePoint]->z;

        Eigen::Vector3d RR(x, y, z);

        element3_d[eno]->givenNormal = RR.cross(RL);
        ;

        double normalMagnitude = element3_d[eno]->givenNormal.norm();
        element3_d[eno]->area = 0.5 * normalMagnitude;
        element3_d[eno]->givenNormal /= normalMagnitude;

        // access to element vertex
        x1 = gvertex[element3_d[eno]->vertexIds[0]]->x;
        y1 = gvertex[element3_d[eno]->vertexIds[0]]->y;
        z1 = gvertex[element3_d[eno]->vertexIds[0]]->z;

        x2 = gvertex[element3_d[eno]->vertexIds[1]]->x;
        y2 = gvertex[element3_d[eno]->vertexIds[1]]->y;
        z2 = gvertex[element3_d[eno]->vertexIds[1]]->z;

        x3 = gvertex[element3_d[eno]->vertexIds[2]]->x;
        y3 = gvertex[element3_d[eno]->vertexIds[2]]->y;
        z3 = gvertex[element3_d[eno]->vertexIds[2]]->z;

        // calculate element mid-point
        element3_d[eno]->midPoint(0, 0) = (x1 + x2 + x3) / 3.0;
        element3_d[eno]->midPoint(1, 0) = (y1 + y2 + y3) / 3.0;
        element3_d[eno]->midPoint(2, 0) = (z1 + z2 + z3) / 3.0;

        // inertia utilities
        x_2 = ((x1 + x2) * (x2 + x3) + x1 * x1 + x3 * x3) / 12.0;
        y_2 = ((y1 + y2) * (y2 + y3) + y1 * y1 + y3 * y3) / 12.0;
        z_2 = ((z1 + z2) * (z2 + z3) + z1 * z1 + z3 * z3) / 12.0;

        xy = ((x1 + x2 + x3) * (y1 + y2 + y3) + x1 * y1 + x2 * y2 + x3 * y3) / 24.0;
        xz = ((x1 + x2 + x3) * (z1 + z2 + z3) + x1 * z1 + x2 * z2 + x3 * z3) / 24.0;
        yz = ((y1 + y2 + y3) * (z1 + z2 + z3) + y1 * z1 + y2 * z2 + y3 * z3) / 24.0;

        x_3 = ((x1 + x2 + x3) * (x1 * x1 + x2 * x2 + x3 * x3) + x1 * x2 * x3) / 20.0;
        y_3 = ((y1 + y2 + y3) * (y1 * y1 + y2 * y2 + y3 * y3) + y1 * y2 * y3) / 20.0;
        z_3 = ((z1 + z2 + z3) * (z1 * z1 + z2 * z2 + z3 * z3) + z1 * z2 * z3) / 20.0;

        x_2y = ((3.0 * y1 + y2 + y3) * x1 * x1 + (y1 + 3.0 * y2 + y3) * x2 * x2 + (y1 + y2 + 3.0 * y3) * x3 * x3
                + (2.0 * y1 + 2.0 * y2 + y3) * x1 * x2 + (2.0 * y1 + y2 + 2.0 * y3) * x1 * x3
                + (y1 + 2.0 * y2 + 2.0 * y3) * x2 * x3)
               / 60.0;
        x_2z = ((3.0 * z1 + z2 + z3) * x1 * x1 + (z1 + 3.0 * z2 + z3) * x2 * x2 + (z1 + z2 + 3.0 * z3) * x3 * x3
                + (2.0 * z1 + 2.0 * z2 + z3) * x1 * x2 + (2.0 * z1 + z2 + 2.0 * z3) * x1 * x3
                + (z1 + 2.0 * z2 + 2.0 * z3) * x2 * x3)
               / 60.0;

        y_2x = ((3.0 * x1 + x2 + x3) * y1 * y1 + (x1 + 3.0 * x2 + x3) * y2 * y2 + (x1 + x2 + 3.0 * x3) * y3 * y3
                + (2.0 * x1 + 2.0 * x2 + x3) * y1 * y2 + (2.0 * x1 + x2 + 2.0 * x3) * y1 * y3
                + (x1 + 2.0 * x2 + 2.0 * x3) * y2 * y3)
               / 60.0;
        y_2z = ((3.0 * z1 + z2 + z3) * y1 * y1 + (z1 + 3.0 * z2 + z3) * y2 * y2 + (z1 + z2 + 3.0 * z3) * y3 * y3
                + (2.0 * z1 + 2.0 * z2 + z3) * y1 * y2 + (2.0 * z1 + z2 + 2.0 * z3) * y1 * y3
                + (z1 + 2.0 * z2 + 2.0 * z3) * y2 * y3)
               / 60.0;

        z_2y = ((3.0 * y1 + y2 + y3) * z1 * z1 + (y1 + 3.0 * y2 + y3) * z2 * z2 + (y1 + y2 + 3.0 * y3) * z3 * z3
                + (2.0 * y1 + 2.0 * y2 + y3) * z1 * z2 + (2.0 * y1 + y2 + 2.0 * y3) * z1 * z3
                + (y1 + 2.0 * y2 + 2.0 * y3) * z2 * z3)
               / 60.0;
        z_2x = ((3.0 * x1 + x2 + x3) * z1 * z1 + (x1 + 3.0 * x2 + x3) * z2 * z2 + (x1 + x2 + 3.0 * x3) * z3 * z3
                + (2.0 * x1 + 2.0 * x2 + x3) * z1 * z2 + (2.0 * x1 + x2 + 2.0 * x3) * z1 * z3
                + (x1 + 2.0 * x2 + 2.0 * x3) * z2 * z3)
               / 60.0;

        xyz = ((x1 + x2 + x3) * (y1 + y2 + y3) * (z1 + z2 + z3) - (y2 * z3 + y3 * z2 - 4.0 * y1 * z1) * x1 / 2.0
               - (y1 * z3 + y3 * z1 - 4.0 * y2 * z2) * x2 / 2.0 - (y1 * z2 + y2 * z1 - 4.0 * y3 * z3) * x3 / 2.0)
              / 60.0;

        xRow1.push_back(Eigen::Vector3d(x_2, 2 * xy, 2 * xz));
        yRow1.push_back(Eigen::Vector3d(2 * xy, y_2, 2 * yz));
        zRow1.push_back(Eigen::Vector3d(2 * xz, 2 * yz, z_2));

        xRow2.push_back(Eigen::Vector3d(x_2y, y_2x, 2 * xyz));
        yRow2.push_back(Eigen::Vector3d(x_2z, 2 * xyz, z_2x));
        zRow2.push_back(Eigen::Vector3d(2 * xyz, y_2z, z_2y));

        xRow3.push_back(Eigen::Vector3d(x_3, 3 * x_2y, 3 * x_2z));
        yRow3.push_back(Eigen::Vector3d(3 * y_2x, y_3, 3 * y_2z));
        zRow3.push_back(Eigen::Vector3d(3 * z_2x, 3 * z_2y, z_3));
    }

    double m000(0.0), m100(0.0), m010(0.0), m001(0.0);
    double m110(0.0), m101(0.0), m011(0.0), m200(0.0), m020(0.0), m002(0.0);

    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        // zeroth order moment (same as volume enclosed by the mesh)
        m000 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(element3_d[eno]->midPoint);

        // first order moments (together they specify the centroid of the mesh)
        m100 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(xRow1[eno]);
        m010 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(yRow1[eno]);
        m001 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(zRow1[eno]);

        m110 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(xRow2[eno]);
        m101 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(yRow2[eno]);
        m011 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(zRow2[eno]);

        m200 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(xRow3[eno]);
        m020 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(yRow3[eno]);
        m002 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(zRow3[eno]);
    }

    m000 /= 3.0;

    m100 /= 3.0;
    m010 /= 3.0;
    m001 /= 3.0;

    m110 /= 3.0;
    m101 /= 3.0;
    m011 /= 3.0;

    m200 /= 4.5;
    m020 /= 4.5;
    m002 /= 4.5;

    // calculate the inertia tensor
    double Ixx, Iyy, Izz, Ixy, Ixz, Iyz;
    Ixx = (m020 + m002) - ((m010 * m010 + m001 * m001) / m000);
    Iyy = (m200 + m002) - ((m100 * m100 + m001 * m001) / m000);
    Izz = (m200 + m020) - ((m100 * m100 + m010 * m010) / m000);
    Ixy = m110 - m100 * m010 / m000;
    Ixz = m101 - m100 * m001 / m000;
    Iyz = m011 - m010 * m001 / m000;

    // area-weighted moment of inertia
    cm << m100 / m000, m010 / m000, m001 / m000;

    Eigen::Matrix3d matrix{ { Ixx, -Ixy, -Ixz }, { -Ixy, Iyy, -Iyz }, { -Ixz, -Iyz, Izz } };
    Eigen::Matrix3d evecMat{ Eigen::Matrix3d::Zero() };
    Eigen::Vector3d eVal{ Eigen::Vector3d::Zero() };

    eigen_decomposition(matrix, evecMat, eVal);

    Eigen::Vector3d radii;
    for (size_t i = 0; i < 3; ++i)
    {
        radii(i, 0) = 1.0 / sqrt(fabs(eVal[i]));
    }

    // sort principle axis
    vector<float> orig;
    for (size_t i = 0; i < 3; ++i)
    {
        orig.push_back(radii(i, 0));
    }

    // this is a vector of {value,index} pairs
    vector<pair<float, size_t>> vp;
    vp.reserve(orig.size());

    for (size_t i = 0; i != orig.size(); ++i)
    {
        vp.push_back(make_pair(orig[i], i));
    }

    // sorting will put lower values ahead of larger ones, resolving ties using the original index
    sort(vp.begin(), vp.end());

    for (size_t i = 0; i < 3; ++i)
    {
        pcx(i, 0) = evecMat(i, vp[0].second);
        pbx(i, 0) = evecMat(i, vp[1].second);
        pax(i, 0) = evecMat(i, vp[2].second);
    }

    // change the global coordiantes w.r.t. centroid
    for (size_t gno = 0; gno < gvertex.size(); ++gno)
    {
        gvertex[gno]->x -= cm(0, 0);
        gvertex[gno]->y -= cm(1, 0);
        gvertex[gno]->z -= cm(2, 0);
    }

    // make the surface normals (to each triangle elements) oriented outwards
    double dsum(0.0);
    for (size_t i = 0; i < 3; ++i)
    {
        x = gvertex[element3_d[0]->side[i].orientation[1]]->x - gvertex[element3_d[0]->side[i].orientation[0]]->x;
        y = gvertex[element3_d[0]->side[i].orientation[1]]->y + gvertex[element3_d[0]->side[i].orientation[0]]->y;
        dsum += x * y;
    }

    // sign of the area of the triangle is negative (<0), i.e. triangle vertex orientation clock-wise, normal needs to
    // be flipped
    if (dsum > 0)
    {
        for (size_t eno = 0; eno < element3_d.size(); eno++)
        {
            element3_d[eno]->givenNormal *= -1.0;
        }
    }
}

void edgeDescriptors(vector<Vertices*>& gvertex, vector<ElementList*>& e_list, vector<Triangle3D*>& element3_d)
{
    // edge descriptors

    int eno, neno, p1, p2;
    double x, y, z, signal, edgeAngle;

    // calculate the 3D-edge angles and angle-axis pairs
    for (size_t i = 0; i < e_list.size(); ++i)
    {
        x = gvertex[e_list[i]->nodes[0]]->x - gvertex[e_list[i]->nodes[1]]->x;
        y = gvertex[e_list[i]->nodes[0]]->y - gvertex[e_list[i]->nodes[1]]->y;
        z = gvertex[e_list[i]->nodes[0]]->z - gvertex[e_list[i]->nodes[1]]->z;

        // length of an edge
        e_list[i]->length = sqrt(x * x + y * y + z * z);

        eno = e_list[i]->sharedElement[0];
        neno = e_list[i]->sharedElement[1];

        p2 = element3_d[eno]->side[e_list[i]->sharedEdge[0]].excludePoint;
        p1 = element3_d[neno]->side[e_list[i]->sharedEdge[1]].excludePoint;

        x = gvertex[p1]->x - gvertex[p2]->x;
        y = gvertex[p1]->y - gvertex[p2]->y;
        z = gvertex[p1]->z - gvertex[p2]->z;

        Eigen::Vector3d e(x, y, z);
        e /= e.norm();

        // this demands the triangle normal MUST BE OUTWARD !!!
        signal = element3_d[eno]->givenNormal.dot(e);

        double acosCorr = element3_d[eno]->givenNormal.dot(element3_d[neno]->givenNormal);
        edgeAngle = acos(min(max(acosCorr, -1.0), 1.0));

        // assign edge bending angle
        element3_d[eno]->side[e_list[i]->sharedEdge[0]].edgBendAngle = edgeAngle;
        element3_d[neno]->side[e_list[i]->sharedEdge[1]].edgBendAngle = edgeAngle;

        // edge angle convex
        if (signal < 0)
        {
            edgeAngle = PI - edgeAngle;
        }

        // edge angle concave
        else
        {
            edgeAngle = PI + edgeAngle;
        }

        e_list[i]->edge3DAngle = edgeAngle;

        // assign triangle edge-3d-angle (PI +/- edge-bending-angle)
        element3_d[eno]->side[e_list[i]->sharedEdge[0]].edgAngle = e_list[i]->edge3DAngle;
        element3_d[neno]->side[e_list[i]->sharedEdge[1]].edgAngle = e_list[i]->edge3DAngle;
    }

    /// Calculate perimeter for each triangle elements
    for (size_t eno = 0; eno < element3_d.size(); ++eno)
    {
        double l(0.0);

        for (size_t sid = 0; sid < 3; sid++)
        {
            l += e_list[element3_d[eno]->sideToEdge[sid]]->length;

            double sidePcat(element3_d[eno]->side[sid].edgBendAngle / element3_d[eno]->edgAngNorm);

            element3_d[eno]->side[sid].pCat = sidePcat;
        }
        element3_d[eno]->periMeter = l;
    }
}

void connectWithGlobe(vector<Vertices*>& gvertex, vector<Triangle3D*>& element3_d)
{
    // link the global coordinates to the triangles
    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        element3_d[eno]->Vertex[0] << gvertex[element3_d[eno]->vertexIds[0]]->x, gvertex[element3_d[eno]->vertexIds[0]]->y,
        gvertex[element3_d[eno]->vertexIds[0]]->z;
        element3_d[eno]->Vertex[1] << gvertex[element3_d[eno]->vertexIds[1]]->x, gvertex[element3_d[eno]->vertexIds[1]]->y,
        gvertex[element3_d[eno]->vertexIds[1]]->z;
        element3_d[eno]->Vertex[2] << gvertex[element3_d[eno]->vertexIds[2]]->x, gvertex[element3_d[eno]->vertexIds[2]]->y,
        gvertex[element3_d[eno]->vertexIds[2]]->z;

        // mid-point coordinate for each triangle
        element3_d[eno]->midPoint
        = (element3_d[eno]->Vertex[0] + element3_d[eno]->Vertex[1] + element3_d[eno]->Vertex[2]) / 3.0;
    }
}

void Image3dTo2D(vector<Region*>& regions, vector<Triangle3D*>& element3_d, double p_cat_curv_list[])
{
    // rotate all the 3d-triangles and put them on x-y plane (normal z-axis)
    Eigen::Vector3d z;
    z << 0.0, 0.0, 1.0;

    Eigen::Vector3d cross;
    double cosAngle(0.0);
    double crossMagnitude(0.0);

    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        // quarternion rotation angle (cosine)
        cosAngle = element3_d[eno]->givenNormal.dot(z);

        // Eigen::Quaternion rotation axis
        cross = element3_d[eno]->givenNormal.cross(z);
        crossMagnitude = cross.norm();

        Eigen::Vector3d nz(0.0, 0.0, 0.0);

        // assign quarternion rotation angle to this triangle region
        regions[eno]->rotAngle = acos(cosAngle);

        // the 3d-triangle required a rotation (triangle normal is non-parallel to z-axis)
        if ((crossMagnitude > 0) || (cosAngle < 0))
        {
            // triangle normal is not anti-parallel to z-axis
            if (crossMagnitude > 0)
            {
                // nomalize quarternion rotation axis
                cross /= crossMagnitude;
                cross *= sin(regions[eno]->rotAngle / 2);
            }

            // triangle normal is anti-parallel to z-axis
            else
            {
                cross << 0.0, -1.0, 0.0;
            }

            // get the quarternion of this triangle region
            regions[eno]->Q.w() = cos(regions[eno]->rotAngle / 2);
            regions[eno]->Q.x() = cross(0, 0);
            regions[eno]->Q.y() = cross(1, 0);
            regions[eno]->Q.z() = cross(2, 0);

            // find the necessary Eigen::Quaternion for 3d -> 2d triangle rotation and rotate the 3d-triangle to
            // 2d-triangle
            Eigen::Quaternion<double> qr(0.0, 0.0, 0.0, 0.0);

            // rotate the mid-point of the triangle region
            Eigen::Quaternion<double> qMp(
            0, element3_d[eno]->midPoint(0, 0), element3_d[eno]->midPoint(1, 0), element3_d[eno]->midPoint(2, 0));
            qr = regions[eno]->Q * qMp * regions[eno]->Q.inverse();
            regions[eno]->zOffset = qr.z();

            // rotate all other three (=3) vertices of the triangle region
            Eigen::Quaternion<double> qA(
            0, element3_d[eno]->Vertex[0](0, 0), element3_d[eno]->Vertex[0](1, 0), element3_d[eno]->Vertex[0](2, 0));
            qr = regions[eno]->Q * qA * regions[eno]->Q.inverse();
            regions[eno]->vertices[0](0, 0) = qr.x();
            regions[eno]->vertices[0](1, 0) = qr.y();

            Eigen::Quaternion<double> qB(
            0, element3_d[eno]->Vertex[1](0, 0), element3_d[eno]->Vertex[1](1, 0), element3_d[eno]->Vertex[1](2, 0));
            qr = regions[eno]->Q * qB * regions[eno]->Q.inverse();
            regions[eno]->vertices[1](0, 0) = qr.x();
            regions[eno]->vertices[1](1, 0) = qr.y();

            Eigen::Quaternion<double> qC(
            0, element3_d[eno]->Vertex[2](0, 0), element3_d[eno]->Vertex[2](1, 0), element3_d[eno]->Vertex[2](2, 0));
            qr = regions[eno]->Q * qC * regions[eno]->Q.inverse();
            regions[eno]->vertices[2](0, 0) = qr.x();
            regions[eno]->vertices[2](1, 0) = qr.y();

            // rotate the triangle normal it self (to check the rotation accuracy: should be parallel to z-axis )
            Eigen::Quaternion<double> nC(
            0, element3_d[eno]->givenNormal(0, 0), element3_d[eno]->givenNormal(1, 0), element3_d[eno]->givenNormal(2, 0));
            qr = regions[eno]->Q * nC * regions[eno]->Q.inverse();
            nz << qr.x(), qr.y(), qr.z();
        }

        // there is no need to rotate the 3d-triangle (triangle normal is parallel to z-axis)
        else
        {
            regions[eno]->Q.w() = 0.0;
            regions[eno]->Q.x() = 0.0;
            regions[eno]->Q.y() = 0.0;
            regions[eno]->Q.z() = 0.0;

            regions[eno]->vertices[0](0, 0) = element3_d[eno]->Vertex[0](0, 0);
            regions[eno]->vertices[0](1, 0) = element3_d[eno]->Vertex[0](1, 0);
            regions[eno]->vertices[1](0, 0) = element3_d[eno]->Vertex[1](0, 0);
            regions[eno]->vertices[1](1, 0) = element3_d[eno]->Vertex[1](1, 0);
            regions[eno]->vertices[2](0, 0) = element3_d[eno]->Vertex[2](0, 0);
            regions[eno]->vertices[2](1, 0) = element3_d[eno]->Vertex[2](1, 0);

            regions[eno]->zOffset = element3_d[eno]->midPoint(2, 0);
            nz = element3_d[eno]->givenNormal;
        }

        // check the accuracy of the Eigen::Quaternion rotation
        nz /= nz.norm();
        double s = nz.dot(z);

        if (s <= 0)
        {
            cout << "Normals not aligned towards z (anti-parallel to) axis !!! emergency exit" << endl;
            exit(-1);
        }

        // change the local coordinate frame to the Center of Mass frame
        regions[eno]->midPoint
        = (regions[eno]->vertices[0] + regions[eno]->vertices[1] + regions[eno]->vertices[2]) / 3.0;
        for (size_t i = 0; i < 3; i++)
        {
            regions[eno]->vertices[i] -= regions[eno]->midPoint;
        }

        // assign the necessary components from 3d-triangle to 2d-triangle
        regions[eno]->regionId = element3_d[eno]->elementId;
        regions[eno]->side = element3_d[eno]->side;
        regions[eno]->sideMap = element3_d[eno]->sideMap;
        regions[eno]->sideRevMap = element3_d[eno]->sideRevMap;
        regions[eno]->area = element3_d[eno]->area;
        regions[eno]->periMeter = element3_d[eno]->periMeter;
        regions[eno]->faceTag = element3_d[eno]->faceTag;

        // assign appropriate pCat at the edges of the triangular region
        double totalPcat(element3_d[eno]->pcatEdgeWeight * p_cat_curv_list[element3_d[eno]->edgeTag]);

        // keep all the probabilities conserved (i.e in the range: 0.0 - 1.0)
        if (totalPcat > 1.0)
        {
            totalPcat = 1.0;
        }

        for (int sid = 0; sid < 3; sid++)
        {
            regions[eno]->side[sid].pCat = totalPcat * element3_d[eno]->side[sid].pCat;
        }
    }
    Eigen::Matrix2d A1, A2;
    Eigen::Vector2d b1, b2;
    int v0, v1, v2;
    double cosTheta, sinTheta, x, y;

    // calculate the rotation angle between two 2d-triangle side
    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            int ceno = regions[eno]->sideMap[j];
            int cesd = regions[ceno]->sideRevMap[regions[eno]->regionId];

            v1 = regions[eno]->side[j].dir[0];
            v2 = regions[eno]->side[j].dir[1];
            regions[eno]->side[j].midPoint = (regions[eno]->vertices[v1] + regions[eno]->vertices[v2]) / 2.0;

            v0 = regions[eno]->side[j].excludePointMap[regions[eno]->side[j].excludePoint];

            A1 << (regions[eno]->vertices[v1](0, 0) - regions[eno]->vertices[v0](0, 0)),
            (regions[eno]->vertices[v2](0, 0) - regions[eno]->vertices[v0](0, 0)),
            (regions[eno]->vertices[v1](1, 0) - regions[eno]->vertices[v0](1, 0)),
            (regions[eno]->vertices[v2](1, 0) - regions[eno]->vertices[v0](1, 0));

            b1 << regions[eno]->vertices[v0](0, 0), regions[eno]->vertices[v0](1, 0);

            Eigen::Vector2d AB = regions[eno]->vertices[v1] - regions[eno]->vertices[v2];
            AB /= AB.norm();

            regions[eno]->side[j].dir2D = AB;
            v1 = regions[ceno]->side[cesd].dir[1];
            v2 = regions[ceno]->side[cesd].dir[0];

            v0 = regions[ceno]->side[cesd].excludePointMap[regions[ceno]->side[cesd].excludePoint];

            A2 << (regions[ceno]->vertices[v1](0, 0) - regions[ceno]->vertices[v0](0, 0)),
            (regions[ceno]->vertices[v2](0, 0) - regions[ceno]->vertices[v0](0, 0)),
            (regions[ceno]->vertices[v1](1, 0) - regions[ceno]->vertices[v0](1, 0)),
            (regions[ceno]->vertices[v2](1, 0) - regions[ceno]->vertices[v0](1, 0));

            b2 << regions[ceno]->vertices[v0](0, 0), regions[ceno]->vertices[v0](1, 0);

            regions[eno]->side[j].A = A2 * A1.inverse();
            regions[eno]->side[j].b = b2 - regions[eno]->side[j].A * b1;

            regions[ceno]->side[cesd].A = A1 * A2.inverse();
            regions[ceno]->side[cesd].b = b1 - regions[ceno]->side[cesd].A * b2;

            Eigen::MatrixXd G = A1.inverse();

            if (fabs(A1.determinant() * A2.determinant()) <= 0.0)
            {
                cout << "Affine matrix with Determinant = 0.0, is detected !!! emergency exit" << endl;
                exit(-1);
            }

            Eigen::Vector2d BA = regions[ceno]->vertices[v1] - regions[ceno]->vertices[v2];
            BA /= BA.norm();

            Eigen::Vector3d s1(AB(0, 0), AB(1, 0), 0.0), s2(BA(0, 0), BA(1, 0), 0.0), ts(0.0, 0.0, 0.0);

            ts = s1.cross(s2);

            cosTheta = s1.dot(s2);
            sinTheta = ts.norm();

            Eigen::Vector3d z(0.0, 0.0, 1.0);
            if (ts.norm() > 0)
            {
                int k;
                ts /= ts.norm();
                k = (int)round(z.dot(ts));
                sinTheta *= k;
            }

            double angle = atan2(sinTheta, cosTheta);

            if (angle < 0)
            {
                angle += 2 * PI;
            }
            if (angle > 2 * PI)
            {
                angle -= 2 * PI;
            }

            regions[eno]->side[j].xyRotationAngle = angle;
        }

        // get the orientation axis (required for calculating the order parameter)
        vector<Eigen::Vector3d> ends3D;
        for (int p = 0; p < 2; p++)
        {
            ends3D.push_back(Eigen::Vector3d(0, 0, 0));
        }

        for (int i = 0; i < 3; i++)
        {
            regions[eno]->orientation.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
        }

        x = 1.0;
        y = 0.0;

        // first two (=2) base vactor : rotate 2-d (x-axis and y-axis) -> triangle normal plane
        for (size_t ax = 0; ax < 2; ax++)
        {
            cross << regions[eno]->Q.x(), regions[eno]->Q.y(), regions[eno]->Q.z();
            crossMagnitude = cross.norm();

            // if rotation required
            if (crossMagnitude > 0.0)
            {
                Eigen::Quaternion<double> qr(0.0, 0.0, 0.0, 0.0);

                // rotate the base vector from 2d-plane to 3d-triangle plane
                Eigen::Quaternion<double> qI0(0, 0.0, 0.0, 0.0);
                qr = regions[eno]->Q.inverse() * qI0 * regions[eno]->Q;
                ends3D[0](0, 0) = qr.x();
                ends3D[0](1, 0) = qr.y();
                ends3D[0](2, 0) = qr.z();

                Eigen::Quaternion<double> qI1(0, x, y, 0.0);
                qr = regions[eno]->Q.inverse() * qI1 * regions[eno]->Q;
                ends3D[1](0, 0) = qr.x();
                ends3D[1](1, 0) = qr.y();
                ends3D[1](2, 0) = qr.z();
            }

            // if no rotation required
            else
            {
                ends3D[0](0, 0) = 0.0;
                ends3D[0](1, 0) = 0.0;
                ends3D[0](2, 0) = 0.0;

                ends3D[1](0, 0) = x;
                ends3D[1](1, 0) = y;
                ends3D[1](2, 0) = 0.0;
            }

            regions[eno]->orientation[ax] << (ends3D[1](0, 0) - ends3D[0](0, 0)), (ends3D[1](1, 0) - ends3D[0](1, 0)),
            (ends3D[1](2, 0) - ends3D[0](2, 0));

            // change base vector from x-axis[(0,0),(1,0)] to y-axis[(0,0),(0,1)]
            swap(x, y);
        }

        // 3rd base vector (z-axis -> triangle normal map)
        regions[eno]->orientation[2] = element3_d[eno]->givenNormal;

        ends3D.clear();
    }
}

void viewGraph(vector<Vertices*>& gvertex, vector<Triangle3D*>& element3_d, Eigen::Vector3d& centroid, string output_dir)
{
    string filename1(output_dir + "/ViewMeshDefault.off");
    string filename2(output_dir + "/ViewMeshFaceTag.dat");
    string filename3(output_dir + "/ViewMeshEdgeTag.dat");
    fstream fp1, fp2, fp3;

    fp1.open(filename1.c_str(), ios::out);
    fp2.open(filename2.c_str(), ios::out);
    fp3.open(filename3.c_str(), ios::out);

    fp1 << "OFF" << endl;
    fp1 << gvertex.size() << " " << element3_d.size() << " " << 0 << endl;

    for (size_t i = 0; i < gvertex.size(); i++)
    {
        fp1 << gvertex[i]->x + centroid(0, 0) << " " << gvertex[i]->y + centroid(1, 0) << " "
            << gvertex[i]->z + centroid(2, 0) << endl;
    }

    for (size_t i = 0; i < element3_d.size(); i++)
    {
        fp1 << 3 << " ";
        for (size_t j = 0; j < 3; j++)
        {
            fp1 << element3_d[i]->vertexIds[j] << " ";
        }
        fp1 << endl;

        fp2 << element3_d[i]->faceTag << " ";
        fp3 << element3_d[i]->edgeTag << " ";
    }

    fp1.close();
    fp2.close();
    fp3.close();
}

bool linePlaneIntersect(Eigen::Vector3d n, Eigen::Vector3d p, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d& i)
{
    Eigen::Vector3d ba = b - a;
    float nDotA = n.dot(a);
    float nDotBA = n.dot(ba);
    float d = n.dot(p);

    float t = (d - nDotA) / nDotBA;
    bool intersect(false);

    if (t >= 0.0f && t <= 1.0f)
    {
        intersect = true;
        i = a + t * ba;
    }

    return (intersect);
}

double areaPolygon3D(vector<Vertices>& v, Eigen::Vector3d n)
{
    size_t v_size(v.size());
    v.push_back(Vertices(v[0].x, v[0].y, v[0].z));

    float area = 0;
    float an, ax, ay, az;
    int coord;
    size_t i, j, k;

    if (v_size < 3)
    {
        return 0;
    }

    // select largest abs coordinate to ignore for projection
    ax = (n(0, 0) > 0 ? n(0, 0) : -n(0, 0));
    ay = (n(1, 0) > 0 ? n(1, 0) : -n(1, 0));
    az = (n(2, 0) > 0 ? n(2, 0) : -n(2, 0));

    coord = 3;
    if (ax > ay)
    {
        if (ax > az)
        {
            coord = 1;
        }
    }
    else if (ay > az)
    {
        coord = 2;
    }

    // compute area of the 2D projection
    switch (coord)
    {
        case 1:
            for (i = 1, j = 2, k = 0; i < v_size; i++, j++, k++)
            {
                area += (v[i].y * (v[j].z - v[k].z));
            }
            break;

        case 2:
            for (i = 1, j = 2, k = 0; i < v_size; i++, j++, k++)
            {
                area += (v[i].z * (v[j].x - v[k].x));
            }
            break;

        case 3:
            for (i = 1, j = 2, k = 0; i < v_size; i++, j++, k++)
            {
                area += (v[i].x * (v[j].y - v[k].y));
            }
            break;
    }

    switch (coord)
    {
        case 1:
            area += (v[v_size].y * (v[1].z - v[v_size - 1].z));
            break;
        case 2:
            area += (v[v_size].z * (v[1].x - v[v_size - 1].x));
            break;
        case 3:
            area += (v[v_size].x * (v[1].y - v[v_size - 1].y));
            break;
    }

    // scale to get area before projection
    an = sqrt(ax * ax + ay * ay + az * az);
    switch (coord)
    {
        case 1:
            area *= (an / (2 * n(0, 0)));
            break;
        case 2:
            area *= (an / (2 * n(1, 0)));
            break;
        case 3:
            area *= (an / (2 * n(2, 0)));
    }

    return fabs(area);
}

double intersectingPolygon(vector<Vertices>& polygon,
                           vector<Vertices>& gvertex,
                           vector<Region*>& regions,
                           vector<ElementList>& e_list,
                           Eigen::Vector3d np,
                           Eigen::Vector3d cp,
                           string shape_type)
{
    double polyArea(0.0);
    if (shape_type == "2D-plane")
    {
        cout << "2d-plane, so no division plane" << endl;
    }
    else
    {
        polygon.clear();

        for (size_t eno = 0; eno < regions.size(); eno++)
        {
            // FOLLOW-UP
            // This used to be a comparison (regions[eno]->polyIntersectMark==0).
            // Check if the fix introduces any regressions.
            regions[eno]->polyIntersectMark = 0;
        }

        vector<int> intersectedEdges;
        double x, y, z;
        for (size_t i = 0; i < e_list.size(); i++)
        {
            x = gvertex[e_list[i].nodes[0]].x;
            y = gvertex[e_list[i].nodes[0]].y;
            z = gvertex[e_list[i].nodes[0]].z;
            Eigen::Vector3d lS(x, y, z);

            x = gvertex[e_list[i].nodes[1]].x;
            y = gvertex[e_list[i].nodes[1]].y;
            z = gvertex[e_list[i].nodes[1]].z;
            Eigen::Vector3d lE(x, y, z);

            Eigen::Vector3d i_vec(0.0, 0.0, 0.0);

            bool Intersect = linePlaneIntersect(np, cp, lS, lE, i_vec);

            if (Intersect)
            {
                intersectedEdges.push_back(i);
                e_list[i].intersectionByPlane = i_vec;
            }
        }

        for (size_t i = 0; i < intersectedEdges.size(); i++)
        {
            int edg = intersectedEdges[i];

            int faceID = e_list[edg].sharedElement[0];
            regions[faceID]->intersectEdg.push_back(edg);
            regions[faceID]->polyIntersectMark = 1;

            faceID = e_list[edg].sharedElement[1];
            regions[faceID]->intersectEdg.push_back(edg);
            regions[faceID]->polyIntersectMark = 1;
        }

        vector<int> polygonEdg;
        int startEdgeIndex = intersectedEdges[0];
        int currentEdgeIndex = startEdgeIndex;

        polygonEdg.push_back(currentEdgeIndex);
        int currentFace = e_list[currentEdgeIndex].sharedElement[0];

        while (true)
        {
            currentEdgeIndex = (currentEdgeIndex == regions[currentFace]->intersectEdg[0])
                               ? regions[currentFace]->intersectEdg[1]
                               : regions[currentFace]->intersectEdg[0];

            currentFace = (currentFace == e_list[currentEdgeIndex].sharedElement[0])
                          ? e_list[currentEdgeIndex].sharedElement[1]
                          : e_list[currentEdgeIndex].sharedElement[0];

            if (currentEdgeIndex == startEdgeIndex)
            {
                break;
            }

            polygonEdg.push_back(currentEdgeIndex);
        }

        // construct polygon
        for (size_t i = 0; i < polygonEdg.size(); i++)
        {
            int edg = polygonEdg[i];
            Eigen::Vector3d pv = e_list[edg].intersectionByPlane;
            polygon.push_back(Vertices(pv(0, 0), pv(1, 0), pv(2, 0)));
        }

        // polygon area
        polyArea = areaPolygon3D(polygon, np);

        intersectedEdges.clear();
        polygonEdg.clear();

        for (size_t i = 0; i < regions.size(); i++)
        {
            regions[i]->intersectEdg.clear();
        }
    }
    return (polyArea);
}

vector<string> split(string str, char delimiter)
{
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while (getline(ss, tok, delimiter))
    {
        internal.push_back(tok);
    }

    return internal;
}

/// PBC for 2D-plane and disc
void establishPBC(vector<Vertices*>& vertices,
                  vector<Triangle3D*>& triangles,
                  vector<Triangle3D*>& view_triangles,
                  vector<ElementList*>& e_list,
                  double size_factor)
{
    int bCount(0);
    vector<int> edgeInMesh, edgeInBoundary;

    /// find the boundary edges
    for (size_t i = 0; i < e_list.size(); i++)
    {
        bool hitorigin(e_list[i]->nodes[1] == 5);

        if (hitorigin)
        {
            bCount++;
        }

        if (!hitorigin)
        {
            edgeInMesh.push_back(i);
        }
    }

    int firstBT = triangles.size() - bCount;

    for (size_t i = 0; i < edgeInMesh.size(); i++)
    {
        int n1 = e_list[edgeInMesh[i]]->sharedElement[0];
        int n2 = e_list[edgeInMesh[i]]->sharedElement[1];
        int largest = ((n1 > n2) ? n1 : n2);

        bool inBoundary(largest >= firstBT);

        if (inBoundary)
        {
            edgeInBoundary.push_back(edgeInMesh[i]);
        }
    }

    edgeInMesh.clear();

    /// delete PPB communicating triangles
    for (int i = 0; i < bCount; i++)
    {
        triangles.pop_back();
        view_triangles.pop_back();
    }

    /// boundary edge (clock/anti-clock wise) sorting
    int vid(vertices.size() - 1);
    for (size_t i = 0; i < edgeInBoundary.size(); i++)
    {
        vid++;
        e_list[edgeInBoundary[i]]->midPoint = vid;

        double xm
        = (vertices[e_list[edgeInBoundary[i]]->nodes[0]]->x + vertices[e_list[edgeInBoundary[i]]->nodes[1]]->x) / 2.0;
        double ym
        = (vertices[e_list[edgeInBoundary[i]]->nodes[0]]->y + vertices[e_list[edgeInBoundary[i]]->nodes[1]]->y) / 2.0;
        double zm
        = (vertices[e_list[edgeInBoundary[i]]->nodes[0]]->z + vertices[e_list[edgeInBoundary[i]]->nodes[1]]->z) / 2.0;

        Eigen::Vector3d vm(xm, ym, zm);
        vertices.push_back(new Vertices(vm(0, 0), vm(1, 0), vm(2, 0)));
    }

    vector<double> edgePositionBoundary;

    for (size_t i = 0; i < edgeInBoundary.size(); i++)
    {
        Eigen::Vector2d radialExt(vertices[e_list[edgeInBoundary[i]]->midPoint]->x,
                                  vertices[e_list[edgeInBoundary[i]]->midPoint]->y);
        radialExt /= radialExt.norm();

        double theta = atan2(radialExt(1, 0), radialExt(0, 0));

        if (theta < 0)
        {
            theta += 2 * PI;
        }

        if (theta > 2 * PI)
        {
            theta -= 2 * PI;
        }

        edgePositionBoundary.push_back(theta);
    }

    /// remove un-necessary vertices
    for (size_t i = 0; i <= edgeInBoundary.size(); i++)
    {
        vertices.pop_back();
    }

    vector<pair<double, int>> vp;
    vp.reserve(edgePositionBoundary.size());

    for (size_t i = 0; i != edgePositionBoundary.size(); i++)
    {
        vp.push_back(make_pair(edgePositionBoundary[i], i));
    }
    edgePositionBoundary.clear();

    sort(vp.begin(), vp.end());

    vector<int> indexBoundary;
    for (size_t i = 0; i < vp.size(); i++)
    {
        indexBoundary.push_back(vp[i].second);
    }

    int n1, n2, p1, p2, e1, e2, s1, s2, l1, l2;

    /// establish top-botton PBC connectivity
    e1 = edgeInBoundary[indexBoundary[0]];
    n1 = e_list[e1]->sharedElement[0];
    n2 = e_list[e1]->sharedElement[1];
    s1 = e_list[e1]->sharedEdge[0];
    s2 = e_list[e1]->sharedEdge[1];
    if (n1 < n2)
    {
        p1 = n1;
        l1 = s1;
    }
    else
    {
        p1 = n2;
        l1 = s2;
    }
    e2 = edgeInBoundary[indexBoundary[2]];
    n1 = e_list[e2]->sharedElement[0];
    n2 = e_list[e2]->sharedElement[1];
    s1 = e_list[e2]->sharedEdge[0];
    s2 = e_list[e2]->sharedEdge[1];
    if (n1 < n2)
    {
        p2 = n1;
        l2 = s1;
    }
    else
    {
        p2 = n2;
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
    e1 = edgeInBoundary[indexBoundary[1]];
    n1 = e_list[e1]->sharedElement[0];
    n2 = e_list[e1]->sharedElement[1];
    s1 = e_list[e1]->sharedEdge[0];
    s2 = e_list[e1]->sharedEdge[1];
    if (n1 < n2)
    {
        p1 = n1;
        l1 = s1;
    }
    else
    {
        p1 = n2;
        l1 = s2;
    }
    e2 = edgeInBoundary[indexBoundary[3]];
    n1 = e_list[e2]->sharedElement[0];
    n2 = e_list[e2]->sharedElement[1];
    s1 = e_list[e2]->sharedEdge[0];
    s2 = e_list[e2]->sharedEdge[1];
    if (n1 < n2)
    {
        p2 = n1;
        l2 = s1;
    }
    else
    {
        p2 = n2;
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
