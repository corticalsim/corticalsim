#include "mesh.hpp"
#include "geometry.hpp"
#include "region.hpp"
#include "system.hpp"
#include "mt_tip.hpp"
#include "linalg.hpp"

// initialize edge
void ini_edge_record(struct EdgeRecord* head, vector<int>& v, int cid, int fid)
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
void add_edge_record(struct EdgeRecord* head, vector<int>& v, int cid, int fid)
{
    EdgeRecord* newedge_record = new EdgeRecord;
    newedge_record->key = v;
    newedge_record->cell_id = cid;
    newedge_record->edge_index = fid;
    newedge_record->next = nullptr;
    EdgeRecord* cur = head;

    while (cur)
    {
        if (cur->next == nullptr)
        {
            cur->next = newedge_record;
            return;
        }
        cur = cur->next;
    }
}

// delete from edge front
void delete_first(struct EdgeRecord** head)
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
    int list_count(0);
    while (list)
    {
        cout << "-----------------------------" << endl;
        cout << "Member: " << list_count << endl;
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
        list_count++;
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
        EdgeRecord** bin;
        bin = new EdgeRecord*[node_max];

        for (int j = 0; j < node_max; j++)
        {
            bin[j] = new EdgeRecord;
            ini_edge_record(bin[j], v1, -1, -1);
        }

        EdgeRecord* edge = *e_list;
        while (edge)
        {
            add_edge_record(bin[edge->key[i]], edge->key, edge->cell_id, edge->edge_index);

            EdgeRecord* oldedge = edge;
            edge = edge->next;
            delete oldedge;
            oldedge = nullptr;
        }

        for (int j = 0; j < node_max; j++)
        {
            delete_first(&bin[j]);
        }

        *e_list = new EdgeRecord;
        ini_edge_record(*e_list, v1, -1, -1);

        for (int j = 0; j < node_max; j++)
        {
            EdgeRecord* binedge = bin[j];
            while (binedge)
            {
                add_edge_record(*e_list, binedge->key, binedge->cell_id, binedge->edge_index);
                binedge = binedge->next;
            }
        }
        delete_first(&(*e_list));
        delete[] bin;
        bin = nullptr;
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
    string cell_index(g->system->p.geomParam);

    // get show mesh condition
    int show_mesh(g->system->p.showMesh);

    // get name of the input triangle mesh to be imported
    string in_file(g->system->p.geometry);

    // scale the mesh, if asked
    double size_factor(g->system->p.geometryScaleFactor);

    // global vertex list
    vector<Vertices*> gvertex;

    // all triangle elemet list
    vector<Triangle3D*> element3_d;

    // centroid/principle-axis of the cell
    Eigen::Vector3d centroid(0.0, 0.0, 0.0);
    Eigen::Vector3d principle_axis(0.0, 0.0, 0.0);
    Eigen::Vector3d ax_b(0.0, 0.0, 0.0);
    Eigen::Vector3d ax_c(0.0, 0.0, 0.0);

    // list of triangle-edges to be detected
    vector<ElementList*> e_list;

    int a0(0), a1(0), a2(0);
    int vmax(0), trimax(0), spf(0), cur(0);
    string str;

    string filename;
    double x(0), y(0), z(0), e_norm, e_ang(0.0);

    if (in_file == "cubeReal")
    {
        gvertex.push_back(new Vertices(0.5 * size_factor, 0.5 * size_factor, 1.0 * size_factor)); // 0
        gvertex.push_back(new Vertices(0.5 * size_factor, 0.5 * size_factor, 0.0));              // 1

        gvertex.push_back(new Vertices(0.5 * size_factor, 0.0, 0.5 * size_factor));              // 2
        gvertex.push_back(new Vertices(0.5 * size_factor, 1.0 * size_factor, 0.5 * size_factor)); // 3

        gvertex.push_back(new Vertices(1.0 * size_factor, 0.5 * size_factor, 0.5 * size_factor)); // 4
        gvertex.push_back(new Vertices(0.0, 0.5 * size_factor, 0.5 * size_factor));              // 5

        gvertex.push_back(new Vertices(0.0, 0.0, 1.0 * size_factor));                           // 6
        gvertex.push_back(new Vertices(1.0 * size_factor, 0.0, 1.0 * size_factor));              // 7
        gvertex.push_back(new Vertices(1.0 * size_factor, 1.0 * size_factor, 1.0 * size_factor)); // 8
        gvertex.push_back(new Vertices(0.0, 1.0 * size_factor, 1.0 * size_factor));              // 9

        gvertex.push_back(new Vertices(0.0, 0.0, 0.0));                           // 10
        gvertex.push_back(new Vertices(1.0 * size_factor, 0.0, 0.0));              // 11
        gvertex.push_back(new Vertices(1.0 * size_factor, 1.0 * size_factor, 0.0)); // 12
        gvertex.push_back(new Vertices(0.0, 1.0 * size_factor, 0.0));              // 13

        // top-face
        /// 1st
        element3_d.push_back(new Triangle3D());
        element3_d[0]->vertexIds[0] = 0;
        element3_d[0]->vertexIds[1] = 6;
        element3_d[0]->vertexIds[2] = 7;
        element3_d[0]->faceTag = 1;

        /// 2nd
        element3_d.push_back(new Triangle3D());
        element3_d[1]->vertexIds[0] = 0;
        element3_d[1]->vertexIds[1] = 7;
        element3_d[1]->vertexIds[2] = 8;
        element3_d[1]->faceTag = 1;

        /// 3rd
        element3_d.push_back(new Triangle3D());
        element3_d[2]->vertexIds[0] = 0;
        element3_d[2]->vertexIds[1] = 8;
        element3_d[2]->vertexIds[2] = 9;
        element3_d[2]->faceTag = 1;

        /// 4th
        element3_d.push_back(new Triangle3D());
        element3_d[3]->vertexIds[0] = 0;
        element3_d[3]->vertexIds[1] = 9;
        element3_d[3]->vertexIds[2] = 6;
        element3_d[3]->faceTag = 1;

        // bottom-face
        /// 1st
        element3_d.push_back(new Triangle3D());
        element3_d[4]->vertexIds[0] = 1;
        element3_d[4]->vertexIds[1] = 10;
        element3_d[4]->vertexIds[2] = 11;
        element3_d[4]->faceTag = 1;

        /// 2nd
        element3_d.push_back(new Triangle3D());
        element3_d[5]->vertexIds[0] = 1;
        element3_d[5]->vertexIds[1] = 11;
        element3_d[5]->vertexIds[2] = 12;
        element3_d[5]->faceTag = 1;

        /// 3rd
        element3_d.push_back(new Triangle3D());
        element3_d[6]->vertexIds[0] = 1;
        element3_d[6]->vertexIds[1] = 12;
        element3_d[6]->vertexIds[2] = 13;
        element3_d[6]->faceTag = 1;

        /// 4th
        element3_d.push_back(new Triangle3D());
        element3_d[7]->vertexIds[0] = 1;
        element3_d[7]->vertexIds[1] = 13;
        element3_d[7]->vertexIds[2] = 10;
        element3_d[7]->faceTag = 1;

        // left-face
        /// 1st
        element3_d.push_back(new Triangle3D());
        element3_d[8]->vertexIds[0] = 2;
        element3_d[8]->vertexIds[1] = 6;
        element3_d[8]->vertexIds[2] = 7;
        element3_d[8]->faceTag = 3;

        /// 2nd
        element3_d.push_back(new Triangle3D());
        element3_d[9]->vertexIds[0] = 2;
        element3_d[9]->vertexIds[1] = 7;
        element3_d[9]->vertexIds[2] = 11;
        element3_d[9]->faceTag = 3;

        /// 3rd
        element3_d.push_back(new Triangle3D());
        element3_d[10]->vertexIds[0] = 2;
        element3_d[10]->vertexIds[1] = 11;
        element3_d[10]->vertexIds[2] = 10;
        element3_d[10]->faceTag = 3;

        /// 4th
        element3_d.push_back(new Triangle3D());
        element3_d[11]->vertexIds[0] = 2;
        element3_d[11]->vertexIds[1] = 10;
        element3_d[11]->vertexIds[2] = 6;
        element3_d[11]->faceTag = 3;

        // right-face
        /// 1st
        element3_d.push_back(new Triangle3D());
        element3_d[12]->vertexIds[0] = 3;
        element3_d[12]->vertexIds[1] = 8;
        element3_d[12]->vertexIds[2] = 9;
        element3_d[12]->faceTag = 3;

        /// 2nd
        element3_d.push_back(new Triangle3D());
        element3_d[13]->vertexIds[0] = 3;
        element3_d[13]->vertexIds[1] = 9;
        element3_d[13]->vertexIds[2] = 13;
        element3_d[13]->faceTag = 3;

        /// 3rd
        element3_d.push_back(new Triangle3D());
        element3_d[14]->vertexIds[0] = 3;
        element3_d[14]->vertexIds[1] = 13;
        element3_d[14]->vertexIds[2] = 12;
        element3_d[14]->faceTag = 3;

        /// 4th
        element3_d.push_back(new Triangle3D());
        element3_d[15]->vertexIds[0] = 3;
        element3_d[15]->vertexIds[1] = 12;
        element3_d[15]->vertexIds[2] = 8;
        element3_d[15]->faceTag = 3;

        // front-face
        /// 1st
        element3_d.push_back(new Triangle3D());
        element3_d[16]->vertexIds[0] = 4;
        element3_d[16]->vertexIds[1] = 7;
        element3_d[16]->vertexIds[2] = 8;
        element3_d[16]->faceTag = 4;

        /// 2nd
        element3_d.push_back(new Triangle3D());
        element3_d[17]->vertexIds[0] = 4;
        element3_d[17]->vertexIds[1] = 8;
        element3_d[17]->vertexIds[2] = 12;
        element3_d[17]->faceTag = 4;

        /// 3rd
        element3_d.push_back(new Triangle3D());
        element3_d[18]->vertexIds[0] = 4;
        element3_d[18]->vertexIds[1] = 12;
        element3_d[18]->vertexIds[2] = 11;
        element3_d[18]->faceTag = 4;

        /// 4th
        element3_d.push_back(new Triangle3D());
        element3_d[19]->vertexIds[0] = 4;
        element3_d[19]->vertexIds[1] = 11;
        element3_d[19]->vertexIds[2] = 7;
        element3_d[19]->faceTag = 4;

        // back-face
        /// 1st
        element3_d.push_back(new Triangle3D());
        element3_d[20]->vertexIds[0] = 5;
        element3_d[20]->vertexIds[1] = 6;
        element3_d[20]->vertexIds[2] = 9;
        element3_d[20]->faceTag = 4;

        /// 2nd
        element3_d.push_back(new Triangle3D());
        element3_d[21]->vertexIds[0] = 5;
        element3_d[21]->vertexIds[1] = 9;
        element3_d[21]->vertexIds[2] = 13;
        element3_d[21]->faceTag = 4;

        /// 3rd
        element3_d.push_back(new Triangle3D());
        element3_d[22]->vertexIds[0] = 5;
        element3_d[22]->vertexIds[1] = 13;
        element3_d[22]->vertexIds[2] = 10;
        element3_d[22]->faceTag = 4;

        /// 4th
        element3_d.push_back(new Triangle3D());
        element3_d[23]->vertexIds[0] = 5;
        element3_d[23]->vertexIds[1] = 10;
        element3_d[23]->vertexIds[2] = 6;
        element3_d[23]->faceTag = 4;

        for (size_t i = 0; i < element3_d.size(); i++)
        {
            element3_d[i]->elementId = i;
        }
    }

    else if (in_file == "2D-plane")
    {
        gvertex.push_back(new Vertices(0.0, 0.0, 0.0));
        gvertex.push_back(new Vertices(0.0, 1.0 * size_factor, 0.0));
        gvertex.push_back(new Vertices(1.0 * size_factor, 1.0 * size_factor, 0.0));
        gvertex.push_back(new Vertices(1.0 * size_factor, 0.0, 0.0));
        gvertex.push_back(new Vertices(0.5 * size_factor, 0.5 * size_factor, 0.0));
        gvertex.push_back(new Vertices(0.5 * size_factor, 0.5 * size_factor, -1.0));

        element3_d.push_back(new Triangle3D());
        element3_d[0]->vertexIds[0] = 0;
        element3_d[0]->vertexIds[1] = 1;
        element3_d[0]->vertexIds[2] = 4;

        element3_d.push_back(new Triangle3D());
        element3_d[1]->vertexIds[0] = 1;
        element3_d[1]->vertexIds[1] = 2;
        element3_d[1]->vertexIds[2] = 4;

        element3_d.push_back(new Triangle3D());
        element3_d[2]->vertexIds[0] = 2;
        element3_d[2]->vertexIds[1] = 3;
        element3_d[2]->vertexIds[2] = 4;

        element3_d.push_back(new Triangle3D());
        element3_d[3]->vertexIds[0] = 3;
        element3_d[3]->vertexIds[1] = 0;
        element3_d[3]->vertexIds[2] = 4;

        element3_d.push_back(new Triangle3D());
        element3_d[4]->vertexIds[0] = 0;
        element3_d[4]->vertexIds[1] = 1;
        element3_d[4]->vertexIds[2] = 5;

        element3_d.push_back(new Triangle3D());
        element3_d[5]->vertexIds[0] = 1;
        element3_d[5]->vertexIds[1] = 2;
        element3_d[5]->vertexIds[2] = 5;

        element3_d.push_back(new Triangle3D());
        element3_d[6]->vertexIds[0] = 2;
        element3_d[6]->vertexIds[1] = 3;
        element3_d[6]->vertexIds[2] = 5;

        element3_d.push_back(new Triangle3D());
        element3_d[7]->vertexIds[0] = 3;
        element3_d[7]->vertexIds[1] = 0;
        element3_d[7]->vertexIds[2] = 5;
    }

    else
    {
        // get input mest file path
        string input_mesh_file = g->system->p.inputDir + "/" + g->system->p.geometry + "/" + in_file + ".dat";

        // open input mesh file
        fstream fp;
        fp.open(input_mesh_file.c_str(), ios::in);

        fp >> vmax >> trimax;

        // read the global vertex coordinates
        for (int i = 0; i < vmax; i++)
        {
            fp >> x >> y >> z;
            gvertex.push_back(new Vertices(0.0, 0.0, 0.0));

            // scale the mesh, if asked
            gvertex[i]->x = x * size_factor;
            gvertex[i]->y = y * size_factor;
            gvertex[i]->z = z * size_factor;
        }

        // read triangle informatio:
        // (a) vertex coordinates
        // (b) spcial face tag
        // (c) curvature tag
        // (d) Pcat norm factor
        // (d) average edge angle (of the curvature domain it belongs)
        for (int i = 0; i < trimax; i++)
        {
            fp >> a0 >> a1 >> a2 >> spf >> cur >> e_norm >> e_ang;
            element3_d.push_back(new Triangle3D());

            element3_d[i]->elementId = i;
            element3_d[i]->faceTag = spf;
            element3_d[i]->edgeTag = cur;
            element3_d[i]->edgAngNorm = e_norm;
            element3_d[i]->pcatEdgeWeight = e_ang;

            // vertices need zero-based numbering
            element3_d[i]->vertexIds[0] = a0 - 1;
            element3_d[i]->vertexIds[1] = a1 - 1;
            element3_d[i]->vertexIds[2] = a2 - 1;
        }

        fp >> str;

        // close the input mesh file
        fp.close();
    }

    // establish connectivity among the triangles
    make_graph(gvertex.size(), e_list, element3_d);

    // store the numer of triangles, vertices and edges
    int f = element3_d.size();
    int v = gvertex.size();
    int e = e_list.size();

    // check mesh consistency
    bool convex(check_eular(f, v, e));

    // copy triangle elements for viewing mesh (visualization in meshLab)
    vector<Triangle3D*> view_triangle;
    for (size_t i = 0; i < element3_d.size(); i++)
    {
        view_triangle.push_back(new Triangle3D());
        view_triangle[i]->vertexIds[0] = element3_d[i]->vertexIds[0];
        view_triangle[i]->vertexIds[1] = element3_d[i]->vertexIds[1];
        view_triangle[i]->vertexIds[2] = element3_d[i]->vertexIds[2];

        view_triangle[i]->faceTag = element3_d[i]->faceTag;
        view_triangle[i]->edgeTag = element3_d[i]->edgeTag;
    }

    // check surface orientation
    orient_surface(gvertex, element3_d, view_triangle);

    // check surface orientation
    bool surface_oriented(check_surface_orientation(element3_d));

    // measure rigidbody properties
    rigid_body_properties(gvertex, element3_d, centroid, principle_axis, ax_b, ax_c);

    // assign coordinate values to the triangle vertices
    connect_with_globe(gvertex, element3_d);

    // get triangle edge related measures/properties
    edge_descriptors(gvertex, e_list, element3_d);

    /// PBC for 2D-plane
    if (in_file == "2D-plane")
    {
        principle_axis << 1.0, 0.0, 0.0;
        ax_b << 0.0, 1.0, 0.0;
        ax_c << 0.0, 0.0, 0.0;
        establish_pbc(gvertex, element3_d, view_triangle, e_list, size_factor);
    }

    // visualize mesh in meshLab
    view_graph(gvertex, view_triangle, centroid, g->system->p.outputDir + "/" + g->system->p.geometry + "/outputData");

    double total_area(0.0);
    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        total_area += element3_d[eno]->area;
    }

    // show mesh in enabled
    if (show_mesh == 1)
    {
        cout << "============================================" << endl;

        // show principle axis
        cout << "Principle Axis ---> (" << principle_axis(0, 0) << "," << principle_axis(1, 0) << ","
             << principle_axis(2, 0) << ")" << endl;

        // show the numer of triangles, vertices and edges
        cout << "F = " << f << " V = " << v << " E = " << e << endl;

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

        if (surface_oriented)
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
    fp << "TriangleNumbers\t" << element3_d.size() << endl;
    fp << "Centroid\t" << centroid(0, 0) << "," << centroid(1, 0) << "," << centroid(2, 0) << endl;
    fp << "PA\t" << principle_axis(0, 0) << "," << principle_axis(1, 0) << "," << principle_axis(2, 0) << endl;
    fp << "PB\t" << ax_b(0, 0) << "," << ax_b(1, 0) << "," << ax_b(2, 0) << endl;
    fp << "PC\t" << ax_c(0, 0) << "," << ax_c(1, 0) << "," << ax_c(2, 0) << endl;
    fp << "Area\t" << total_area << endl;

    fp.close();

    // starting simulation on the surface-mesh
    map<string, double>::iterator it;

    // list of edge catastrophe values
    int edg_number = g->system->p.edgCatMap.size();
    double p_cat_list[edg_number];
    for (int i = 0; i <= edg_number; i++)
    {
        p_cat_list[i] = 0.0;
    }

    for (it = g->system->p.edgCatMap.begin(); it != g->system->p.edgCatMap.end(); it++)
    {
        vector<string> curr_edg = split(it->first, '_');
        int edg_idx{ std::stoi(curr_edg[1]) };
        double edg_cat_val = it->second;
        p_cat_list[edg_idx] = edg_cat_val;
    }

    // list of face catastrophe values
    int face_number = g->system->p.faceCatMap.size();
    for (int i = 0; i <= face_number; i++)
    {
        g->system->p.RegionKcatMultiplier.push_back(1.0);
    }

    for (it = g->system->p.faceCatMap.begin(); it != g->system->p.faceCatMap.end(); it++)
    {
        vector<string> curr_face = split(it->first, '_');
        int face_idx{ std::stoi(curr_face[1]) };
        double face_cat_val = it->second;
        g->system->p.RegionKcatMultiplier[face_idx] = face_cat_val;
    }

    // assign the centroid and principle axis to the geometry
    g->objectCM = centroid;
    g->objectPA = principle_axis;
    g->nucleousPosition << 0.0, 0.0, 0.0;

    // total surface area of the geometry
    g->area = total_area;

    // create only a single path area (at first place)
    g->patchArea[0] = g->area;

    // include all triangles to a single path (at first place)
    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        g->RegionsIndex[0].push_back(eno);
    }

    // total number of region in the geometry
    g->elementMax = element3_d.size();

    // to create PPB
    // (a)  global vertex to be used to create PPB
    for (size_t i = 0; i < gvertex.size(); i++)
    {
        g->globalVertex.push_back(Vertices(0.0, 0.0, 0.0));
        g->globalVertex[i].x = gvertex[i]->x;
        g->globalVertex[i].y = gvertex[i]->y;
        g->globalVertex[i].z = gvertex[i]->z;
    }

    // (b) element list
    for (size_t i = 0; i < e_list.size(); i++)
    {
        g->ppbEdgeList.push_back(ElementList());

        for (size_t j = 0; j < 2; j++)
        {
            g->ppbEdgeList[i].nodes[j] = e_list[i]->nodes[j];
            g->ppbEdgeList[i].sharedElement[j] = e_list[i]->sharedElement[j];
            g->ppbEdgeList[i].sharedEdge[j] = e_list[i]->sharedEdge[j];
        }
    }

    // assign memory for all the 2d-triangular regions of the geometry
    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        g->regions.push_back(new Triangle(element3_d[eno]->area, g));
    }

    // copy 3d-triangle properties information to 2d-triangular regions
    image3d_to2_d(g->regions, element3_d, p_cat_list);

    // 3d-triangle properties are transfered to 2d-triangular regions, so release all the memory occupied by 3d
    // triangles
    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        delete element3_d[eno];
        // FOLLOW-UP
        // This used to be a comparison (element3D[eno] == NULL).
        // Check if the fix introduces any regressions.
        element3_d[eno] = nullptr;
    }
    element3_d.clear();

    // global vertices no longer needed, so release the memories occupied by the global vertices
    for (size_t gno = 0; gno < gvertex.size(); gno++)
    {
        delete gvertex[gno];
        gvertex[gno] = nullptr;
    }
    gvertex.clear();

    // edge list is no longer needed, so release the memories occupied by the edges
    for (size_t el = 0; el < e_list.size(); el++)
    {
        delete e_list[el];
        e_list[el] = nullptr;
    }
    e_list.clear();

    return;
}

void make_graph(int gvertex, vector<ElementList*>& e_list, vector<Triangle3D*>& element3_d)
{
    // make the connectivity graph
    void edge_sort(EdgeRecord**, int);
    void find_connectivity(EdgeRecord**, vector<ElementList*>&);

    struct EdgeRecord* edg_list = new EdgeRecord();

    vector<int> v;
    v.push_back(0);
    v.push_back(0);

    ini_edge_record(edg_list, v, -1, -1);

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

            add_edge_record(edg_list, v, i, j);
            v.clear();
        }
        element3_d[i]->elementId = i;
    }
    delete_first(&(edg_list));

    edge_sort(&edg_list, gvertex);
    find_connectivity(&edg_list, e_list);

    for (size_t i = 0; i < e_list.size(); i++)
    {
        element3_d[e_list[i]->sharedElement[0]]->sideMap[e_list[i]->sharedEdge[0]] = e_list[i]->sharedElement[1];
        element3_d[e_list[i]->sharedElement[1]]->sideMap[e_list[i]->sharedEdge[1]] = e_list[i]->sharedElement[0];

        element3_d[e_list[i]->sharedElement[0]]->sideRevMap[e_list[i]->sharedElement[1]] = e_list[i]->sharedEdge[0];
        element3_d[e_list[i]->sharedElement[1]]->sideRevMap[e_list[i]->sharedElement[0]] = e_list[i]->sharedEdge[1];

        element3_d[e_list[i]->sharedElement[0]]->sideToEdge[e_list[i]->sharedEdge[0]] = i;
        element3_d[e_list[i]->sharedElement[1]]->sideToEdge[e_list[i]->sharedEdge[1]] = i;
    }

    delete edg_list;
    edg_list = nullptr;

    return;
}

void orient_surface(vector<Vertices*>& gvertex, vector<Triangle3D*>& element3_d, vector<Triangle3D*>& view_triangle)
{
    // orient surface normals
    void correct_orientation(vector<Vertices*>&, Triangle3D*&, Triangle3D*&, int, int, int&);

    vector<Triangle3D*> triangle_alive;
    triangle_alive.push_back(element3_d[0]);

    vector<Triangle3D*> triangle_dead;

    while (!triangle_alive.empty())
    {
        Triangle3D* next_triangle = triangle_alive.back();
        triangle_alive.pop_back();

        triangle_dead.push_back(next_triangle);

        for (int i = 0; i < 3; i++)
        {
            Triangle3D* neigh_triangle = element3_d[next_triangle->sideMap[i]];

            bool is_present = find(triangle_dead.begin(), triangle_dead.end(), neigh_triangle) != triangle_dead.end();
            if (!is_present)
            {
                int flip(1);
                int n1 = neigh_triangle->sideRevMap[next_triangle->elementId];
                correct_orientation(gvertex, next_triangle, neigh_triangle, i, n1, flip);
                triangle_alive.push_back(neigh_triangle);
                if (flip == -1)
                {
                    view_triangle[next_triangle->sideMap[i]]->vertexIds[1]
                    = element3_d[next_triangle->sideMap[i]]->vertexIds[2];
                    view_triangle[next_triangle->sideMap[i]]->vertexIds[2]
                    = element3_d[next_triangle->sideMap[i]]->vertexIds[1];
                }
            }
        }
    }
}

void correct_orientation(
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

bool check_surface_orientation(vector<Triangle3D*>& element3_d)
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

bool check_eular(int f, int v, int e)
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

void rigid_body_properties(vector<Vertices*>& gvertex,
                         vector<Triangle3D*>& element3_d,
                         Eigen::Vector3d& cm,
                         Eigen::Vector3d& pax,
                         Eigen::Vector3d& pbx,
                         Eigen::Vector3d& pcx)
{
    // rigidBodyProperties
    vector<Eigen::Vector3d> x_row1, y_row1, z_row1, x_row2, y_row2, z_row2, x_row3, y_row3, z_row3;
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
        Eigen::Vector3d rl(x, y, z);

        x = gvertex[element3_d[eno]->side[0].orientation[0]]->x - gvertex[element3_d[eno]->side[0].excludePoint]->x;
        y = gvertex[element3_d[eno]->side[0].orientation[0]]->y - gvertex[element3_d[eno]->side[0].excludePoint]->y;
        z = gvertex[element3_d[eno]->side[0].orientation[0]]->z - gvertex[element3_d[eno]->side[0].excludePoint]->z;

        Eigen::Vector3d rr(x, y, z);

        element3_d[eno]->givenNormal = rr.cross(rl);
        ;

        double normal_magnitude = element3_d[eno]->givenNormal.norm();
        element3_d[eno]->area = 0.5 * normal_magnitude;
        element3_d[eno]->givenNormal /= normal_magnitude;

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

        x_row1.push_back(Eigen::Vector3d(x_2, 2 * xy, 2 * xz));
        y_row1.push_back(Eigen::Vector3d(2 * xy, y_2, 2 * yz));
        z_row1.push_back(Eigen::Vector3d(2 * xz, 2 * yz, z_2));

        x_row2.push_back(Eigen::Vector3d(x_2y, y_2x, 2 * xyz));
        y_row2.push_back(Eigen::Vector3d(x_2z, 2 * xyz, z_2x));
        z_row2.push_back(Eigen::Vector3d(2 * xyz, y_2z, z_2y));

        x_row3.push_back(Eigen::Vector3d(x_3, 3 * x_2y, 3 * x_2z));
        y_row3.push_back(Eigen::Vector3d(3 * y_2x, y_3, 3 * y_2z));
        z_row3.push_back(Eigen::Vector3d(3 * z_2x, 3 * z_2y, z_3));
    }

    double m000(0.0), m100(0.0), m010(0.0), m001(0.0);
    double m110(0.0), m101(0.0), m011(0.0), m200(0.0), m020(0.0), m002(0.0);

    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        // zeroth order moment (same as volume enclosed by the mesh)
        m000 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(element3_d[eno]->midPoint);

        // first order moments (together they specify the centroid of the mesh)
        m100 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(x_row1[eno]);
        m010 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(y_row1[eno]);
        m001 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(z_row1[eno]);

        m110 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(x_row2[eno]);
        m101 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(y_row2[eno]);
        m011 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(z_row2[eno]);

        m200 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(x_row3[eno]);
        m020 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(y_row3[eno]);
        m002 += element3_d[eno]->area * element3_d[eno]->givenNormal.dot(z_row3[eno]);
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
    double ixx, iyy, izz, ixy, ixz, iyz;
    ixx = (m020 + m002) - ((m010 * m010 + m001 * m001) / m000);
    iyy = (m200 + m002) - ((m100 * m100 + m001 * m001) / m000);
    izz = (m200 + m020) - ((m100 * m100 + m010 * m010) / m000);
    ixy = m110 - m100 * m010 / m000;
    ixz = m101 - m100 * m001 / m000;
    iyz = m011 - m010 * m001 / m000;

    // area-weighted moment of inertia
    cm << m100 / m000, m010 / m000, m001 / m000;

    Eigen::Matrix3d matrix{ { ixx, -ixy, -ixz }, { -ixy, iyy, -iyz }, { -ixz, -iyz, izz } };
    Eigen::Matrix3d evec_mat{ Eigen::Matrix3d::Zero() };
    Eigen::Vector3d e_val{ Eigen::Vector3d::Zero() };

    eigen_decomposition(matrix, evec_mat, e_val);

    Eigen::Vector3d radii;
    for (size_t i = 0; i < 3; ++i)
    {
        radii(i, 0) = 1.0 / sqrt(fabs(e_val[i]));
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
        pcx(i, 0) = evec_mat(i, vp[0].second);
        pbx(i, 0) = evec_mat(i, vp[1].second);
        pax(i, 0) = evec_mat(i, vp[2].second);
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

void edge_descriptors(vector<Vertices*>& gvertex, vector<ElementList*>& e_list, vector<Triangle3D*>& element3_d)
{
    // edge descriptors

    int eno, neno, p1, p2;
    double x, y, z, signal, edge_angle;

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

        double acos_corr = element3_d[eno]->givenNormal.dot(element3_d[neno]->givenNormal);
        edge_angle = acos(min(max(acos_corr, -1.0), 1.0));

        // assign edge bending angle
        element3_d[eno]->side[e_list[i]->sharedEdge[0]].edgBendAngle = edge_angle;
        element3_d[neno]->side[e_list[i]->sharedEdge[1]].edgBendAngle = edge_angle;

        // edge angle convex
        if (signal < 0)
        {
            edge_angle = PI - edge_angle;
        }

        // edge angle concave
        else
        {
            edge_angle = PI + edge_angle;
        }

        e_list[i]->edge3DAngle = edge_angle;

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

            double side_pcat(element3_d[eno]->side[sid].edgBendAngle / element3_d[eno]->edgAngNorm);

            element3_d[eno]->side[sid].pCat = side_pcat;
        }
        element3_d[eno]->periMeter = l;
    }
}

void connect_with_globe(vector<Vertices*>& gvertex, vector<Triangle3D*>& element3_d)
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

void image3d_to2_d(vector<Region*>& regions, vector<Triangle3D*>& element3_d, double p_cat_curv_list[])
{
    // rotate all the 3d-triangles and put them on x-y plane (normal z-axis)
    Eigen::Vector3d z;
    z << 0.0, 0.0, 1.0;

    Eigen::Vector3d cross;
    double cos_angle(0.0);
    double cross_magnitude(0.0);

    for (size_t eno = 0; eno < element3_d.size(); eno++)
    {
        // quarternion rotation angle (cosine)
        cos_angle = element3_d[eno]->givenNormal.dot(z);

        // Eigen::Quaternion rotation axis
        cross = element3_d[eno]->givenNormal.cross(z);
        cross_magnitude = cross.norm();

        Eigen::Vector3d nz(0.0, 0.0, 0.0);

        // assign quarternion rotation angle to this triangle region
        regions[eno]->rotAngle = acos(cos_angle);

        // the 3d-triangle required a rotation (triangle normal is non-parallel to z-axis)
        if ((cross_magnitude > 0) || (cos_angle < 0))
        {
            // triangle normal is not anti-parallel to z-axis
            if (cross_magnitude > 0)
            {
                // nomalize quarternion rotation axis
                cross /= cross_magnitude;
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
            Eigen::Quaternion<double> q_mp(
            0, element3_d[eno]->midPoint(0, 0), element3_d[eno]->midPoint(1, 0), element3_d[eno]->midPoint(2, 0));
            qr = regions[eno]->Q * q_mp * regions[eno]->Q.inverse();
            regions[eno]->zOffset = qr.z();

            // rotate all other three (=3) vertices of the triangle region
            Eigen::Quaternion<double> q_a(
            0, element3_d[eno]->Vertex[0](0, 0), element3_d[eno]->Vertex[0](1, 0), element3_d[eno]->Vertex[0](2, 0));
            qr = regions[eno]->Q * q_a * regions[eno]->Q.inverse();
            regions[eno]->vertices[0](0, 0) = qr.x();
            regions[eno]->vertices[0](1, 0) = qr.y();

            Eigen::Quaternion<double> q_b(
            0, element3_d[eno]->Vertex[1](0, 0), element3_d[eno]->Vertex[1](1, 0), element3_d[eno]->Vertex[1](2, 0));
            qr = regions[eno]->Q * q_b * regions[eno]->Q.inverse();
            regions[eno]->vertices[1](0, 0) = qr.x();
            regions[eno]->vertices[1](1, 0) = qr.y();

            Eigen::Quaternion<double> q_c(
            0, element3_d[eno]->Vertex[2](0, 0), element3_d[eno]->Vertex[2](1, 0), element3_d[eno]->Vertex[2](2, 0));
            qr = regions[eno]->Q * q_c * regions[eno]->Q.inverse();
            regions[eno]->vertices[2](0, 0) = qr.x();
            regions[eno]->vertices[2](1, 0) = qr.y();

            // rotate the triangle normal it self (to check the rotation accuracy: should be parallel to z-axis )
            Eigen::Quaternion<double> n_c(
            0, element3_d[eno]->givenNormal(0, 0), element3_d[eno]->givenNormal(1, 0), element3_d[eno]->givenNormal(2, 0));
            qr = regions[eno]->Q * n_c * regions[eno]->Q.inverse();
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
        double total_pcat(element3_d[eno]->pcatEdgeWeight * p_cat_curv_list[element3_d[eno]->edgeTag]);

        // keep all the probabilities conserved (i.e in the range: 0.0 - 1.0)
        if (total_pcat > 1.0)
        {
            total_pcat = 1.0;
        }

        for (int sid = 0; sid < 3; sid++)
        {
            regions[eno]->side[sid].pCat = total_pcat * element3_d[eno]->side[sid].pCat;
        }
    }
    Eigen::Matrix2d a1, a2;
    Eigen::Vector2d b1, b2;
    int v0, v1, v2;
    double cos_theta, sin_theta, x, y;

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

            a1 << (regions[eno]->vertices[v1](0, 0) - regions[eno]->vertices[v0](0, 0)),
            (regions[eno]->vertices[v2](0, 0) - regions[eno]->vertices[v0](0, 0)),
            (regions[eno]->vertices[v1](1, 0) - regions[eno]->vertices[v0](1, 0)),
            (regions[eno]->vertices[v2](1, 0) - regions[eno]->vertices[v0](1, 0));

            b1 << regions[eno]->vertices[v0](0, 0), regions[eno]->vertices[v0](1, 0);

            Eigen::Vector2d ab = regions[eno]->vertices[v1] - regions[eno]->vertices[v2];
            ab /= ab.norm();

            regions[eno]->side[j].dir2D = ab;
            v1 = regions[ceno]->side[cesd].dir[1];
            v2 = regions[ceno]->side[cesd].dir[0];

            v0 = regions[ceno]->side[cesd].excludePointMap[regions[ceno]->side[cesd].excludePoint];

            a2 << (regions[ceno]->vertices[v1](0, 0) - regions[ceno]->vertices[v0](0, 0)),
            (regions[ceno]->vertices[v2](0, 0) - regions[ceno]->vertices[v0](0, 0)),
            (regions[ceno]->vertices[v1](1, 0) - regions[ceno]->vertices[v0](1, 0)),
            (regions[ceno]->vertices[v2](1, 0) - regions[ceno]->vertices[v0](1, 0));

            b2 << regions[ceno]->vertices[v0](0, 0), regions[ceno]->vertices[v0](1, 0);

            regions[eno]->side[j].A = a2 * a1.inverse();
            regions[eno]->side[j].b = b2 - regions[eno]->side[j].A * b1;

            regions[ceno]->side[cesd].A = a1 * a2.inverse();
            regions[ceno]->side[cesd].b = b1 - regions[ceno]->side[cesd].A * b2;

            Eigen::MatrixXd g = a1.inverse();

            if (fabs(a1.determinant() * a2.determinant()) <= 0.0)
            {
                cout << "Affine matrix with Determinant = 0.0, is detected !!! emergency exit" << endl;
                exit(-1);
            }

            Eigen::Vector2d ba = regions[ceno]->vertices[v1] - regions[ceno]->vertices[v2];
            ba /= ba.norm();

            Eigen::Vector3d s1(ab(0, 0), ab(1, 0), 0.0), s2(ba(0, 0), ba(1, 0), 0.0), ts(0.0, 0.0, 0.0);

            ts = s1.cross(s2);

            cos_theta = s1.dot(s2);
            sin_theta = ts.norm();

            Eigen::Vector3d z(0.0, 0.0, 1.0);
            if (ts.norm() > 0)
            {
                int k;
                ts /= ts.norm();
                k = (int)round(z.dot(ts));
                sin_theta *= k;
            }

            double angle = atan2(sin_theta, cos_theta);

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
        vector<Eigen::Vector3d> ends3_d;
        for (int p = 0; p < 2; p++)
        {
            ends3_d.push_back(Eigen::Vector3d(0, 0, 0));
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
            cross_magnitude = cross.norm();

            // if rotation required
            if (cross_magnitude > 0.0)
            {
                Eigen::Quaternion<double> qr(0.0, 0.0, 0.0, 0.0);

                // rotate the base vector from 2d-plane to 3d-triangle plane
                Eigen::Quaternion<double> q_i0(0, 0.0, 0.0, 0.0);
                qr = regions[eno]->Q.inverse() * q_i0 * regions[eno]->Q;
                ends3_d[0](0, 0) = qr.x();
                ends3_d[0](1, 0) = qr.y();
                ends3_d[0](2, 0) = qr.z();

                Eigen::Quaternion<double> q_i1(0, x, y, 0.0);
                qr = regions[eno]->Q.inverse() * q_i1 * regions[eno]->Q;
                ends3_d[1](0, 0) = qr.x();
                ends3_d[1](1, 0) = qr.y();
                ends3_d[1](2, 0) = qr.z();
            }

            // if no rotation required
            else
            {
                ends3_d[0](0, 0) = 0.0;
                ends3_d[0](1, 0) = 0.0;
                ends3_d[0](2, 0) = 0.0;

                ends3_d[1](0, 0) = x;
                ends3_d[1](1, 0) = y;
                ends3_d[1](2, 0) = 0.0;
            }

            regions[eno]->orientation[ax] << (ends3_d[1](0, 0) - ends3_d[0](0, 0)), (ends3_d[1](1, 0) - ends3_d[0](1, 0)),
            (ends3_d[1](2, 0) - ends3_d[0](2, 0));

            // change base vector from x-axis[(0,0),(1,0)] to y-axis[(0,0),(0,1)]
            swap(x, y);
        }

        // 3rd base vector (z-axis -> triangle normal map)
        regions[eno]->orientation[2] = element3_d[eno]->givenNormal;

        ends3_d.clear();
    }
}

void view_graph(vector<Vertices*>& gvertex, vector<Triangle3D*>& element3_d, Eigen::Vector3d& centroid, string output_dir)
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

bool line_plane_intersect(Eigen::Vector3d n, Eigen::Vector3d p, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d& i)
{
    Eigen::Vector3d ba = b - a;
    float n_dot_a = n.dot(a);
    float n_dot_ba = n.dot(ba);
    float d = n.dot(p);

    float t = (d - n_dot_a) / n_dot_ba;
    bool intersect(false);

    if (t >= 0.0f && t <= 1.0f)
    {
        intersect = true;
        i = a + t * ba;
    }

    return (intersect);
}

double area_polygon3_d(vector<Vertices>& v, Eigen::Vector3d n)
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

double intersecting_polygon(vector<Vertices>& polygon,
                           vector<Vertices>& gvertex,
                           vector<Region*>& regions,
                           vector<ElementList>& e_list,
                           Eigen::Vector3d np,
                           Eigen::Vector3d cp,
                           string shape_type)
{
    double poly_area(0.0);
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

        vector<int> intersected_edges;
        double x, y, z;
        for (size_t i = 0; i < e_list.size(); i++)
        {
            x = gvertex[e_list[i].nodes[0]].x;
            y = gvertex[e_list[i].nodes[0]].y;
            z = gvertex[e_list[i].nodes[0]].z;
            Eigen::Vector3d l_s(x, y, z);

            x = gvertex[e_list[i].nodes[1]].x;
            y = gvertex[e_list[i].nodes[1]].y;
            z = gvertex[e_list[i].nodes[1]].z;
            Eigen::Vector3d l_e(x, y, z);

            Eigen::Vector3d i_vec(0.0, 0.0, 0.0);

            bool intersect = line_plane_intersect(np, cp, l_s, l_e, i_vec);

            if (intersect)
            {
                intersected_edges.push_back(i);
                e_list[i].intersectionByPlane = i_vec;
            }
        }

        for (size_t i = 0; i < intersected_edges.size(); i++)
        {
            int edg = intersected_edges[i];

            int face_id = e_list[edg].sharedElement[0];
            regions[face_id]->intersectEdg.push_back(edg);
            regions[face_id]->polyIntersectMark = 1;

            face_id = e_list[edg].sharedElement[1];
            regions[face_id]->intersectEdg.push_back(edg);
            regions[face_id]->polyIntersectMark = 1;
        }

        vector<int> polygon_edg;
        int start_edge_index = intersected_edges[0];
        int current_edge_index = start_edge_index;

        polygon_edg.push_back(current_edge_index);
        int current_face = e_list[current_edge_index].sharedElement[0];

        while (true)
        {
            current_edge_index = (current_edge_index == regions[current_face]->intersectEdg[0])
                               ? regions[current_face]->intersectEdg[1]
                               : regions[current_face]->intersectEdg[0];

            current_face = (current_face == e_list[current_edge_index].sharedElement[0])
                          ? e_list[current_edge_index].sharedElement[1]
                          : e_list[current_edge_index].sharedElement[0];

            if (current_edge_index == start_edge_index)
            {
                break;
            }

            polygon_edg.push_back(current_edge_index);
        }

        // construct polygon
        for (size_t i = 0; i < polygon_edg.size(); i++)
        {
            int edg = polygon_edg[i];
            Eigen::Vector3d pv = e_list[edg].intersectionByPlane;
            polygon.push_back(Vertices(pv(0, 0), pv(1, 0), pv(2, 0)));
        }

        // polygon area
        poly_area = area_polygon3_d(polygon, np);

        intersected_edges.clear();
        polygon_edg.clear();

        for (size_t i = 0; i < regions.size(); i++)
        {
            regions[i]->intersectEdg.clear();
        }
    }
    return (poly_area);
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
void establish_pbc(vector<Vertices*>& vertices,
                  vector<Triangle3D*>& triangles,
                  vector<Triangle3D*>& view_triangles,
                  vector<ElementList*>& e_list,
                  double size_factor)
{
    int b_count(0);
    vector<int> edge_in_mesh, edge_in_boundary;

    /// find the boundary edges
    for (size_t i = 0; i < e_list.size(); i++)
    {
        bool hitorigin(e_list[i]->nodes[1] == 5);

        if (hitorigin)
        {
            b_count++;
        }

        if (!hitorigin)
        {
            edge_in_mesh.push_back(i);
        }
    }

    int first_bt = triangles.size() - b_count;

    for (size_t i = 0; i < edge_in_mesh.size(); i++)
    {
        int n1 = e_list[edge_in_mesh[i]]->sharedElement[0];
        int n2 = e_list[edge_in_mesh[i]]->sharedElement[1];
        int largest = ((n1 > n2) ? n1 : n2);

        bool in_boundary(largest >= first_bt);

        if (in_boundary)
        {
            edge_in_boundary.push_back(edge_in_mesh[i]);
        }
    }

    edge_in_mesh.clear();

    /// delete PPB communicating triangles
    for (int i = 0; i < b_count; i++)
    {
        triangles.pop_back();
        view_triangles.pop_back();
    }

    /// boundary edge (clock/anti-clock wise) sorting
    int vid(vertices.size() - 1);
    for (size_t i = 0; i < edge_in_boundary.size(); i++)
    {
        vid++;
        e_list[edge_in_boundary[i]]->midPoint = vid;

        double xm
        = (vertices[e_list[edge_in_boundary[i]]->nodes[0]]->x + vertices[e_list[edge_in_boundary[i]]->nodes[1]]->x) / 2.0;
        double ym
        = (vertices[e_list[edge_in_boundary[i]]->nodes[0]]->y + vertices[e_list[edge_in_boundary[i]]->nodes[1]]->y) / 2.0;
        double zm
        = (vertices[e_list[edge_in_boundary[i]]->nodes[0]]->z + vertices[e_list[edge_in_boundary[i]]->nodes[1]]->z) / 2.0;

        Eigen::Vector3d vm(xm, ym, zm);
        vertices.push_back(new Vertices(vm(0, 0), vm(1, 0), vm(2, 0)));
    }

    vector<double> edge_position_boundary;

    for (size_t i = 0; i < edge_in_boundary.size(); i++)
    {
        Eigen::Vector2d radial_ext(vertices[e_list[edge_in_boundary[i]]->midPoint]->x,
                                  vertices[e_list[edge_in_boundary[i]]->midPoint]->y);
        radial_ext /= radial_ext.norm();

        double theta = atan2(radial_ext(1, 0), radial_ext(0, 0));

        if (theta < 0)
        {
            theta += 2 * PI;
        }

        if (theta > 2 * PI)
        {
            theta -= 2 * PI;
        }

        edge_position_boundary.push_back(theta);
    }

    /// remove un-necessary vertices
    for (size_t i = 0; i <= edge_in_boundary.size(); i++)
    {
        vertices.pop_back();
    }

    vector<pair<double, int>> vp;
    vp.reserve(edge_position_boundary.size());

    for (size_t i = 0; i != edge_position_boundary.size(); i++)
    {
        vp.push_back(make_pair(edge_position_boundary[i], i));
    }
    edge_position_boundary.clear();

    sort(vp.begin(), vp.end());

    vector<int> index_boundary;
    for (size_t i = 0; i < vp.size(); i++)
    {
        index_boundary.push_back(vp[i].second);
    }

    int n1, n2, p1, p2, e1, e2, s1, s2, l1, l2;

    /// establish top-botton PBC connectivity
    e1 = edge_in_boundary[index_boundary[0]];
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
    e2 = edge_in_boundary[index_boundary[2]];
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
    e1 = edge_in_boundary[index_boundary[1]];
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
    e2 = edge_in_boundary[index_boundary[3]];
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

    index_boundary.clear();
}
