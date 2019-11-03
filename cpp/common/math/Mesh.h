
#ifndef  Mesh_h
#define  Mesh_h

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <string>

//#include "fastmath.h"
//#include "Vec2.h"
//#include "Vec3.h"
//#include "Mat3.h"
//#include "quaternion.h"
#include "raytrace.h"
#include "geom3D.h"


// implementation
#include "integerOps.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

inline uint32_t vertIdjk( uint8_t i, uint8_t j, uint8_t k ){
    uint8_t t;
    if(i>k){t=i; i=k; k=t; }
    if(i>j){t=i; i=j; j=t; }
    if(j>k){t=j; j=k; k=t; }
    return pack32(i,j,k,0);
}

class Mesh{ public:
    int    rendered_shape;
    Vec3d  center; // center of spherical bounding box
    double Rbound; // radius of spherical bounding box
    /*
    int npoints;
    int ntris;
    int npolys;
    int nedges;

    Vec3d  * points;
    Vec3i  * tris;
    int    * polyNs;
    double ** polys;
    Vec2i  * edges;
    */

    std::vector<Vec3d>        points;
    std::vector<Vec3d>        normals;
    std::vector<Vec3i>        triangles;
    std::vector<Polygon*>     polygons;
    //std::vector<double[4]>  planes;
    std::vector<MeshEdge>     edges;

    std::vector<int>        removed_points;
    std::vector<int>        removed_edges;

    Disk3D * disks;  // used for acceleration of raytracing; for each polygon there is one disk

    std::vector<int> tri2poly;

    //std::vector<LineInterval3d> lijsDEBUG;

    // ==== functions

    //int    fromFileOBJ( std::string fname );
    //int    pickVertex ( const Vec3d &ray0, const Vec3d &hRay );
    //double ray        ( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal );

    // =============== W.I.P

    // convex hull
    //void subdivHull  ( const Vec3i& tri, int n, int *selection, Vec3d * from_points );
    //void findHull_oct( int npoints, Vec3d * from_points );

    // ========== Implementations

    /*
    void assign_tri2poly(){
        if( tri2poly ) delete [] tri2poly; tri2poly = new int[triangles.size()];
        for(int itr=0; itr<triangles.size(); itr++){
            for(int ipl=0; ipl<polygons.size(); ipl++){
                Polygon * pl = polygons[ipl];
                for(int i=0; i<triangles.size(); i++){
                    if( pl-> )
                    tri2poly[itr] =
                }
            }
        }
    }
    */

    void   exportEdges ( Vec2i* vts, Vec2i* fcs ){               for(int i=0;i<edges.size();i++){ if(vts)vts[i]=edges[i].verts; if(fcs)fcs[i]=edges[i].faces; }; }
    Vec2i* exportEdgeVs(){ Vec2i* vts = new Vec2i[edges.size()]; for(int i=0;i<edges.size();i++){ vts[i]=edges[i].verts; }; return vts; }
    Vec2i* exportEdgeFs(){ Vec2i* fcs = new Vec2i[edges.size()]; for(int i=0;i<edges.size();i++){ fcs[i]=edges[i].faces; }; return fcs; }

    Vec3d faceCog( int ipl )const {
        Polygon* pl = polygons.at(ipl);
        int n = pl->ipoints.size();
        Vec3d c; c.set(0.0);
        for(int i=0;i<n; i++){ c.add( points[ pl->ipoints[i] ] ); }
        c.mul(1.0d/n);
        return c;
    };

    double polygonArea( int ipl, Vec3d * pl_normal ){
        Polygon* pl = polygons.at(ipl);
        double S = 0.0;
        Vec3d a  = points[ pl->ipoints[0] ];
        Vec3d ab = points[ pl->ipoints[1] ] - a;
        for(int i=2;i<pl->ipoints.size(); i++){
            Vec3d ac  = points[ pl->ipoints[i] ] - a;
            //Vec3d vT  = ac - ab*( ac.dot(ab)/ac.norm() );
            //double Si = vT.dot(ab);
            Vec3d normal; normal.set_cross(ab,ac);
            if(pl_normal) pl_normal->add(normal);
            double Si = 0.5*normal.norm();
            S += Si;
            ab = ac;
        }
        if(pl_normal) pl_normal->mul(0.5/S);
        return S;
    }

    int findEdges( ){
        std::unordered_map<uint64_t,MeshEdge> edge_map;
        for( int i=0; i<polygons.size(); i++ ){
            //if(i>0) break;
            Polygon* pl = polygons[i];
            Vec2i verts;
            int np  = pl->ipoints.size();
            verts.a = pl->ipoints[np-1];
            //printf( "%i %i \n", i, np );
            for(int j=0;j<np;j++){
                verts.b = pl->ipoints[j];
                uint64_t key = symetric_id(verts);
                //uint64_t key = scalar_id(verts);
                int sz0 = edge_map.size();
                MeshEdge& edge = edge_map[key];   // get valid pointer ( if new alocate, if aold take it )
                //printf( "%i :  %i %i %li %i \n", j, verts.a, verts.b, key, edge_map.size()-sz0 );
                if( edge_map.size() > sz0 ){     // new element
                    edge.setVerts(verts.a,verts.b);
                    edge.faces.a = i;
                    edge.faces.b = -1;
                }else{
                    edge.faces.b = i;
                }
                verts.a=verts.b;
            }
        }
        int i=0;
        for(auto kv : edge_map) {
            edges.push_back( kv.second );
            //printf( "%i (%i,%i)(%i,%i)\n", i, kv.second.verts.a, kv.second.verts.b, kv.second.faces.a, kv.second.faces.b );
            i++;
        }
        return i;
    }


    int fromFileOBJ( std::string fname ){
        //std::ofstream infile;
        std::ifstream infile ( fname );
        //myfile << "Writing this to a file.\n";
        std::string line,word;
        //std::stringstream ss;
        if (!infile.is_open()) { std::cout << "cannot open "<<fname<<"\n"; return -1; }
        //int npoints0   = points  .size();
        //int npolygons0 = polygons.size();
        int i=0;
        while ( std::getline(infile,line) ){
            std::cout << line << '\n';
            //ss.str(line);
            std::stringstream ss(line);
            getline(ss, word, ' ');
            if ( word.compare("v")==0 ){
                Vec3d p;
                getline(ss, word, ' '); p.x = stof(word);
                getline(ss, word, ' '); p.y = stof(word);
                getline(ss, word, ' '); p.z = stof(word);
                points.push_back(p);
                i++;
            }else if ( word.compare("f")==0  ){
                Polygon * p = new Polygon();
                polygons.push_back(p);
                while ( getline(ss, word, ' ')) {
                    size_t of;
                    p->ipoints .push_back( stoi(word,&of)-1 );
                    std::cout << word.substr(of) << '\n';
                    p->inormals.push_back( stoi(word.substr(of+2))-1 );
                };
                i++;
            }else if ( word.compare("vn")==0 ){
                Vec3d p;
                getline(ss, word, ' '); p.x = stof(word);
                getline(ss, word, ' '); p.y = stof(word);
                getline(ss, word, ' '); p.z = stof(word);
                normals.push_back(p);
                i++;
            };
        }
        infile.close();
        return i;
    };

    void polygonsToTriangles( bool t2p ){
        //if(t2p){ if( tri2poly ) delete [] tri2poly; tri2poly = new int[triangles.size()]; }
        for(int i=0; i<polygons.size(); i++){
            Polygon * pl = polygons[i];
            Vec3i tri;
            tri.a = pl->ipoints[0];
            tri.b = pl->ipoints[1];
            for(int j=2; j<pl->ipoints.size(); j++){
                tri.c = pl->ipoints[j];
                //printf( "%i %i (%i,%i,%i)\n", i, j, tri.a,tri.b,tri.c );
                triangles.push_back(tri);
                if(t2p) tri2poly.push_back(i);
                tri.b=tri.c;
            }
        }
    }

    void tris2normals( bool normalize ){
        //normals.erase();
        normals.resize(points.size());
        for( Vec3i tri : triangles ){
            Vec3d a,b,c;
            a = points[tri.a];
            b.set_sub(  points[tri.b], a );
            c.set_sub(  points[tri.c], a );
            a.set_cross(b,c);
            a.normalize();
            normals[tri.a].add(a);
            normals[tri.b].add(a);
            normals[tri.c].add(a);
        };
        if(normalize){
            for(int i=0; i<normals.size(); i++){ normals[i].normalize(); }
        }
    }

    int pickVertex( const Vec3d &ray0, const Vec3d &hRay ){
        double r2min=1e+300;
        int imin=0;
        for(int i=0; i<points.size(); i++){
            double t;
            double r2 = rayPointDistance2( ray0, hRay, points[i], t );
            if(r2<r2min){ imin=i; r2min=r2; }
        }
        return imin;
    };

    double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal, int& imin ){
    //int    imin  = 0;
        Vec3d hX,hY;
        hRay.getSomeOrtho( hX, hY );
        double t_min = 1e+300;
        //Vec3d hitpos_min,normal_min;
        imin = -1;
        for(int i=0; i<triangles.size(); i++ ){
            Vec3i itri = triangles[i];
            Vec3d A = points[itri.x];
            Vec3d B = points[itri.y];
            Vec3d C = points[itri.z];
            Vec3d normal_;
            bool inside_;
            //double t = rayTriangle( ray0, hRay, A, B, C, inside_, hitpos_, normal_ );
            double t = rayTriangle2( ray0, hRay, hX, hY, A, B, C, normal_ );
            //printf( "t=%f\n", t );
            inside_ = (t<0.9e+300 )&&(t>0);
            if( inside_ && ( t<t_min ) ){
                t_min  = t;
                normal = normal_;
                imin   = i;
            }
        }
        return t_min;
    };

    /// Topology operations

    void insertEdgeVertex( int ied ){
        int ip = points.size();
        MeshEdge& ed = edges[ied];
        Vec3d p = (points[ed.verts.a] + points[ed.verts.b])*0.5;
        points.push_back(p);
        int ito; Polygon * pl;

        //printf("%i %i\n", ed.verts.a, ed.verts.b );
        // insert to 1
        pl  = polygons[ed.faces.a];
        ito = pl->findEdgeIndex( ed.verts.a, ed.verts.b );
        pl->insertPoint( ip, ito+1 );
        // insert to 2
        pl  = polygons[ed.faces.b];
        ito = pl->findEdgeIndex( ed.verts.a, ed.verts.b );
        pl->insertPoint( ip, ito+1 );
        // edges
        MeshEdge ed2;
        ed2.verts.a = ip;
        ed2.verts.b = ed.verts.b;
        ed .verts.b = ip;
        ed2.faces   = ed.faces;
        edges.push_back(ed2);
    }

    void cleanRemovedPoints( ){
        int reind[ points.size() ];
        int nvalid=0;
        for(int i=0; i<points.size();         i++){ reind[i]=0; };
        for(int i=0; i<removed_points.size(); i++){ reind[removed_points[i]]=-1; }
        for(int i=0; i<points.size();         i++){
            if( reind[i] != -1 ){
                reind[i]      =nvalid;
                points[nvalid]=points[i];
                nvalid++;
            }
        };
        points.resize(nvalid);
        for(MeshEdge& edi: edges ){
            edi.verts.a = reind[ edi.verts.a ];
            edi.verts.b = reind[ edi.verts.b ];
        }
        for(Polygon* pl: polygons ){
            for(int j=0; j<pl->ipoints.size(); j++ ){
                int& ip = pl->ipoints[j];
                ip=reind[ ip ];
            }
        }
    }

    int colapseEdge( int ied ){
        MeshEdge ed = edges[ied];
        removed_points.push_back(ed.verts.b);
        edges.erase(edges.begin()+ied);
        for( MeshEdge& edi: edges ){
            if(edi.verts.a == ed.verts.b) edi.verts.a = ed.verts.a;
            if(edi.verts.b == ed.verts.b) edi.verts.b = ed.verts.a;
        }
        polygons[ed.faces.a]->removePoint(ed.verts.b);
        polygons[ed.faces.b]->removePoint(ed.verts.b);
        for( int i=0; i<polygons.size(); i++ ){ // this operation has cost N^2
            if( (i==ed.faces.a)||(i==ed.faces.b) ) continue;
            Polygon* pl = polygons[i];
            for(int j=0; j<pl->ipoints.size(); j++ ){
                int& ip = pl->ipoints[j];
                if( ip== ed.verts.b) ip=ed.verts.a;
            }
        }
        return ed.verts.a;
    }

    //int bevelVertex( int i ){}

    //int extrudeFace( int i ){}


    //std::unordered_map<std::int,Vec2i>::iterator pair_insert( std::unordered_map<int,Vec2i>& mymap, int key, int val ){
    //    std::unordered_map<std::int,Vec2i>::iterator got = cutted_faces.find( key );

    bool pair_insert( std::unordered_map<int,Vec2i>& mymap, int key, int val ){
        auto got = mymap.find( key );
        if( got == mymap.end() ){
            mymap.insert( { key, (Vec2i){val,-1} } );
            return false;
        }else{
            got->second.b = val;
            return true;
        }
    }

    void fromPlanes(int n, Plane3D* planes ){
        std::unordered_map<uint32_t,int> verts;
        std::vector<int> faceEdges[n];
        int nv =0;
        for(int i=0; i<n; i++){
            //int ie0=edges.size();
            for(int j=i+1; j<n; j++){
                LineInterval3d lij;
                //printf( "===(%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, j, planes[i].normal.x, planes[i].normal.y, planes[i].normal.z,   planes[j].normal.x, planes[j].normal.y, planes[j].normal.z );
                if( !lij.fromPlanes( planes[i].normal, planes[i].iso, planes[j].normal, planes[j].iso ) ){
                    //printf("colinear intersection %i,%i \n", i,j  );
                    continue;
                }
                //lijsDEBUG.push_back(lij);
                MeshEdge eij;
                decltype(verts.find(0)) it;
                uint32_t iv0,iv1;
                eij.verts.a = -1; // DEBUG
                eij.verts.b = -1; // DEBUG
                for( int k=0; k<n; k++ ){
                    if( (k==i)||(k==j) ){ continue; }
                    int side = lij.trim( planes[k].normal, planes[k].iso );
                    if( lij.t0>lij.t1 ){
                        //printf( "   line Vanish (%i,%i,%i) \n", i,j,k );
                        goto lineVanish;
                    };
                    //printf( "side %i \n", side );
                    if     ( side>0 ){ eij.verts.b = k; }
                    else if( side<0 ){ eij.verts.a = k; }
                }
                if( (eij.verts.a<0) or (eij.verts.b<0) ){
                    //printf( "end not trimmed %i %i \n", eij.verts.a, eij.verts.b );
                    continue; // DEBUG
                }
                //lijsDEBUG.push_back(lij);
                // try insert endpoint 1
                eij.faces.a = i;
                eij.faces.b = j;
                //printf( "   edge (%i,%i)(%i,%i) \n", eij.faces.a, eij.faces.b, eij.verts.a, eij.verts.b );
                iv0 = vertIdjk( i, j, eij.verts.a );
                it = verts.find(iv0);
                if(  it == verts.end() ){
                    //nv = points.size();
                    //printf( "   new point (%i,%i,%i) -> %i #%i\n", i, j, eij.verts.a, nv, iv0 );
                    verts.insert({iv0,nv});
                    eij.verts.a = nv;
                    points.push_back( lij.endPoint0() );
                    nv++;
                }else{
                    //printf( "   old point (%i,%i,%i) -> %i #%i\n", i, j, eij.verts.a, it->second, iv0 );
                    eij.verts.a = it->second;
                }
                // try insert endpoint 2
                iv1 = vertIdjk( i, j, eij.verts.b );
                it = verts.find(iv1);
                if(  it == verts.end() ){
                    //nv = points.size();
                    //printf( "   new point (%i,%i,%i) -> %i #%i \n", i, j, eij.verts.b, nv, iv1 );
                    verts.insert({iv1,nv});
                    eij.verts.b = nv;
                    points.push_back( lij.endPoint1() );
                    nv++;
                }else{
                    //printf( "   old point (%i,%i,%i) -> %i #%i\n", i, j, eij.verts.b, it->second, iv1 );
                    eij.verts.b = it->second;
                }
                //eij.verts.a = points.size(); points.push_back( lij.endPoint0() );
                //eij.verts.b = points.size(); points.push_back( lij.endPoint1() );
                faceEdges[i].push_back( edges.size() );
                faceEdges[j].push_back( edges.size() );
                edges.push_back( eij );
            lineVanish:
                ;
            }
            // order edges to make a polygon
            //int ne = edges.size() - ie0;
            
            if( faceEdges[i].size()>0 ){
                //printf( " face: %i nedges %i \n", i, faceEdges[i].size() );
                Polygon * pl = new Polygon();
                //int ieo   = -1;
                int ie    = faceEdges[i][0];
                int iv    = edges[ie].verts.a;
                int ivEnd = edges[ie].verts.b;
                //printf( "%i(%i,%i)\n", ie, iv, ivEnd );
                //for(int iter=0; iter<ne; iter++){
                // todo - we can perhaps sort it by faster algorithms but for small number of edges-per-polygon it is not worth it
                int ndone=0;
                while( iv != ivEnd ){
                    int ivNew = -1;
                    int j;
                    //for( j : faceEdges[i]){
                    for(int jj=0; jj<faceEdges[i].size(); jj++){
                        j = faceEdges[i][jj];
                        if(ie==j) continue;
                        //if((ie==j)||(ie==ieo)) continue;
                        if      (edges[j].verts.a == iv){ ivNew=edges[j].verts.b; break; }
                        else if (edges[j].verts.b == iv){ ivNew=edges[j].verts.a; break; };
                    }
                    if(ivNew>=0){ // anything was found ?
                        pl->ipoints.push_back(iv);
                        pl->iedges .push_back(ie);
                        //ieo = ie;
                        ie  = j;
                        //printf( "ie %i %i iv %i->%i \n", ie, j, iv, ivNew );
                        iv = ivNew;
                    }else                         { printf("error: cannot find %i \n", iv );      return; }
                    if(ndone>=faceEdges[i].size()){ printf("error: cannot close polygon %i \n", iv ); return; }
                    ndone++;
                }
                pl->ipoints.push_back(iv);
                pl->iedges .push_back(ie);
                polygons.push_back(pl);
            }
            /*
            // using map point->edge we can do polygon search in linear time 
            if( faceEdges[i].size()>2 ){
                Vec2i v2e[faceEdges[i].size()];
                for(int ie : faceEdges[i]){
                    edges[ie].vert.a;
                }
                pl->ipoints.push_back(iv);
                pl->iedges .push_back(ie);
                polygons.push_back(pl);
            }
            */
        }
    }

    int cutByPlane( Vec3d dir, double C ){
        int ied=0;
        int ip0 = points.size();
        int ie0 = edges.size();

        std::unordered_map<int,Vec2i> cutted_faces;

        for( MeshEdge& ed : edges ){
            Vec3d& pa = points[ed.verts.a];
            Vec3d& pb = points[ed.verts.b];
            double ta = dir.dot(pa)-C;
            double tb = dir.dot(pb)-C;
            if(  ta*tb < 0 ){
                double t = ta/(ta-tb);
                if(ta<0){
                    t*=-1;
                    removed_points.push_back(ed.verts.a);
                    ed.verts.a = points.size();
                }else{
                    removed_points.push_back(ed.verts.b);
                    ed.verts.b = points.size();
                }
                points.push_back( pa*(1-t) + pb*t );
                pair_insert( cutted_faces, ed.faces.a, ied );
                pair_insert( cutted_faces, ed.faces.b, ied );
                //cutted_faces.insert( ed.faces.a );
                //cutted_faces.insert( ed.faces.b );
            }else{
                if(ta<0){
                    removed_edges.push_back(ied);
                    // TODO - edges between new points
                    // TODO - faces between new points
                    // TODO - Remove also faces
                }
            }
            ied++;
        }
        return ied;
        // TODO - capped cut ?
    }


// =============================
// =============================
// ============================= iMPLEMENTATION
// =============================
// =============================
// WHY THIS FUNCKING IMPLEMENTATION CANNOT BE MOVED TO .cpp ??????
/*    BECAUSE IT SAYS : "||=== Build: test_Mesh in SimpleSimulationEngine (compiler: GNU GCC Compiler) ===|
../../common/math/CMakeFiles/Mesh.dir/Mesh.cpp.o||In function `Mesh::pickVertex(Vec3TYPE<double> const&, Vec3TYPE<double> const&)':|
Mesh.cpp|| multiple definition of `Mesh::pickVertex(Vec3TYPE<double> const&, Vec3TYPE<double> const&)'|
../../common/SDL2OGL/CMakeFiles/SDL2OGL.dir/__/math/Mesh.cpp.o:Mesh.cpp|| first defined here|
../../common/math/CMakeFiles/Mesh.dir/Mesh.cpp.o||In function `Mesh::ray(Vec3TYPE<double> const&, Vec3TYPE<double> const&, Vec3TYPE<double>&)':|
Mesh.cpp|| multiple definition of `Mesh::ray(Vec3TYPE<double> const&, Vec3TYPE<double> const&, Vec3TYPE<double>&)'|
../../common/SDL2OGL/CMakeFiles/SDL2OGL.dir/__/math/Mesh.cpp.o:Mesh.cpp|| first defined here|
../../common/math/CMakeFiles/Mesh.dir/Mesh.cpp.o||In function `Mesh::subdivHull(Vec3TYPE<int> const&, int, int*, Vec3TYPE<double>*)':|
Mesh.cpp|| multiple definition of `Mesh::subdivHull(Vec3TYPE<int> const&, int, int*, Vec3TYPE<double>*)'|
../../common/SDL2OGL/CMakeFiles/SDL2OGL.dir/__/math/Mesh.cpp.o:Mesh.cpp|| first defined here|
../../common/math/CMakeFiles/Mesh.dir/Mesh.cpp.o||In function `Mesh::findHull_oct(int, Vec3TYPE<double>*)':|
Mesh.cpp|| multiple definition of `Mesh::findHull_oct(int, Vec3TYPE<double>*)'|
../../common/SDL2OGL/CMakeFiles/SDL2OGL.dir/__/math/Mesh.cpp.o:Mesh.cpp|| first defined here|
../../common/math/CMakeFiles/Mesh.dir/Mesh.cpp.o||In function `Mesh::fromFileOBJ(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)':|
Mesh.cpp|| multiple definition of `Mesh::fromFileOBJ(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)'|
../../common/SDL2OGL/CMakeFiles/SDL2OGL.dir/__/math/Mesh.cpp.o:Mesh.cpp|| first defined here|
||error: ld returned 1 exit status|
tests/3D/CMakeFiles/test_Mesh.dir/build.make|120|recipe for target 'tests/3D/test_Mesh' failed|
CMakeFiles/Makefile2|3027|recipe for target 'tests/3D/CMakeFiles/test_Mesh.dir/all' failed|
CMakeFiles/Makefile2|3039|recipe for target 'tests/3D/CMakeFiles/test_Mesh.dir/rule' failed|
/home/prokop/git/SimpleSimulationEngine/cpp/Build/tests/3D/Makefile|175|recipe for target 'tests/3D/CMakeFiles/test_Mesh.dir/rule' failed|
||=== Build failed: 15 error(s), 0 warning(s) (0 minute(s), 2 second(s)) ===|
"


WHY ? WHY ? WHY ? WHY ? WHY ?

Maybe because Mesh is compited to  Draw3D object ?????
to do later

*/





//  =======================
//   Convex Hull 3D w.I.P
//  =======================

// http://thomasdiewald.com/blog/?p=1888
// http://mindthenerd.blogspot.cz/2012/05/fastest-convex-hull-algorithm-ever.html

// TO DO :
//  it is not so simple !!!! The edges connecting "tri" vertexes are not guarantied to be part of convex hull
//  if new vertex is found in two neighboring faces, we have to check if the new edge between them is not beyond the original edge
//  => first must develop true mesh edditin functionality

void subdivHull( const Vec3i& tri, int n, int *selection, Vec3d * from_points ){
    Plane3D plane; plane.fromPoints(from_points[tri.a],from_points[tri.b],from_points[tri.c]);
    double distMax=0.0;
    int imax = -1;
    int n_=0;
    int * out_selection = selection + n;
    for(int i=0; i<n; i++){
        int ip = selection[i];
        Vec3d p = points[ip];
        double dist = plane.dist(p);
        //if( dist > 0){
        if( dist > 1e-16 ){  // this way we filter out a,b,c and other degenerate points
            n++;
            out_selection[i] = ip;
            if( dist > distMax ){
                distMax = dist;
                imax = i;
            }
        }
    }
    if ( n_ ){	// last level of recursion ?
        subdivHull( {tri.a,tri.b,imax}, n_, out_selection, from_points );
        subdivHull( {tri.b,tri.c,imax}, n_, out_selection, from_points );
        subdivHull( {tri.c,tri.a,imax}, n_, out_selection, from_points );
    }else{
        //hull_tris.set( tri );
        //ntris++;
        triangles.push_back( tri );
    }
};

void findHull_oct( int npoints, Vec3d * from_points ){
    // find initial octaedra
    int * selection = new int[3*npoints];
    int        imx,ipx,ipy,imy,imz,ipz;
    double     mx, px, py, my, mz, pz;
    for( int i=0; i<npoints; i++ ){
        selection[i] = i;
        Vec3d p = from_points[i];
        if      ( p.x < mx ){ imx=i; mx=p.x; }
        else if ( p.x > px ){ ipx=i; px=p.x; };
        if      ( p.y < my ){ imy=i; my=p.y; }
        else if ( p.y > py ){ ipy=i; py=p.y; };
        if      ( p.z < mz ){ imz=i; mz=p.z; }
        else if ( p.z > pz ){ ipz=i; pz=p.z; };
    }
    subdivHull( {imx,imy,imz},  npoints, selection, from_points );
    subdivHull( {ipx,imy,imz},  npoints, selection, from_points );
    subdivHull( {imx,ipy,imz},  npoints, selection, from_points );
    subdivHull( {ipx,ipy,imz},  npoints, selection, from_points );
    subdivHull( {imx,imy,ipz},  npoints, selection, from_points );
    subdivHull( {ipx,imy,ipz},  npoints, selection, from_points );
    subdivHull( {imx,ipy,ipz},  npoints, selection, from_points );
    subdivHull( {ipx,ipy,ipz},  npoints, selection, from_points );
    delete selection;
};

// Create initial simplex (tetrahedron, 4 vertices). To do this, the 6 Extreme Points [EP], min/max points in X,Y and Z, of the given pointcloud are extracted.
// From those 6 EP the two most distant build the base-line of the base triangle.
// The most distant point of EP to the base line is the 3rd point of the base-triangle.
// To find the pyramids apex, the most distant point to the base-triangle is searched for in the whole point-list.
// Now having 4 points, the inital pyramid can easily be created.

/*
	void findHull_tetra( ){
		// find initial octaedra
		int       imx,ipx;
		double     mx, px;
		for( int i=0; i<n; i++ ){
			Vec3d p = from_points[i];
			if      ( p.x < mx ){ imx=x; mx=p.x; }
			else if ( p.x > px ){ ipx=x; px=p.x; };
		}
		//Vec3d a,b,c;
		//a.set_sub( from_points[imx], from_points[ipx] );
		//a.getSomeOrtho( b, c );
		Vec3d ahat;   ahat. set_sub( from_points[imx], from_points[ipx] ); ahat.normalize();
		Vec3d ahalf;  ahalf.set_add( from_points[imx], from_points[ipx] ); ahalf.mul(0.5);
		for( int i=0; i<n; i++ ){
			Vec3d p;
			p.set_sub( from_points[i], ahalf );
			double cdot = p.dot( ahat );
			p.add_mul( ahat, cdot );
			double r2 = p
			if      ( p.x < mx ){ imx=x; mx=p.x; }
			else if ( p.x > px ){ ipx=x; px=p.x; };
		}

		for( int i=0; i<n; i++ ){
			Vec3d p = from_points[i];
			if      ( p.x < mx ){ imx=x; mx=p.x; }
			else if ( p.x > px ){ ipx=x; px=p.x; };
		}

		subdivHull( {imx,imy,imz},  npoints, selection );
		subdivHull( {ipx,imy,imz},  npoints, selection );
		subdivHull( {imx,ipy,imz},  npoints, selection );
		subdivHull( {ipx,ipy,imz},  npoints, selection );
		subdivHull( {imx,imy,ipz},  npoints, selection );
		subdivHull( {ipx,imy,ipz},  npoints, selection );
		subdivHull( {imx,ipy,ipz},  npoints, selection );
		subdivHull( {ipx,ipy,ipz},  npoints, selection );
	}

	*/

};












#endif



