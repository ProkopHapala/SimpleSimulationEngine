
/**
 * @file Truss.h
 * @brief Truss class stores mechanical structures composed of sticks (girders) and ropes (cables). It can be used for simulating bridges, cranes, spaceships, etc.
 */

#ifndef  Truss_h
#define  Truss_h

#include <vector>
//#include <unordered_map>
#include <cstdio>
#include <cstring>
//#include <string>
#include "raytrace.h"

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"

class GirderParams{  public:
    int n=5;    // number of segments along the girder
    int m=2;    // number of segments 
    int kind=0;
    double widthX = 0.5;
    double widthY = 0.5;
    int kind_long   = 0;  
    int kind_perp   = 1;
    int kind_zigIn  = 2;
    int kind_zigOut = 3;
};


class TrussEdge{ public:
    int a,b;     // indexes of vertices (nodes, points)
    int type;    // type of the edge, e.g. 0 - stick, 1 - rope
    double l0;   // [m] rest length
    double d;    // [m] diameter
    double crossection; // [m^2] crossection area 
    double mass; // total mass of the edge
    double kT;   // stiffness in tension
    double kP;   // stiffness in pressure

    void fromString(char* str){ sscanf ( str, "%i %i %i",&a,&b,&type); a--; b--; }
    void toString  (char* str){ sprintf( str, "%i %i %i", a+1, b+1, type); }
    void print     (         ){  printf(      "%i %i %i", a+1, b+1, type); }
};

class TrussFace{ public:
    int a,b,c;    // indexes of vertices (nodes, points)
    int type;     // type of the edge, e.g. 0 - stick, 1 - rope
    double area;  // area of the face
    double mass;  // total mass of the face

    void fromString(char* str){ sscanf ( str, "%i %i %i",&a,&b,&c,&type); a--; b--; c--; }
    void toString  (char* str){ sprintf( str, "%i %i %i", a+1, b+1, c+1, type); }
    void print     (         ){  printf(      "%i %i %i", a+1, b+1, c+1, type); }
};

class Truss{ public:
    // TODO: should we extend mesh class ?
    int    rendered_shape;
    std::vector<Vec3d>        points;   // what about mass ?
    std::vector<TrussEdge>    edges;    // 
    std::vector<TrussFace>    faces;    // triangles between edges
    //std::vector<Quat4i>     vols;     // volumes between faces
    std::vector<Vec2i>        blocks;
    std::vector<int>          removed_points;
    std::vector<int>          removed_edges;

    inline int addPoint( const Vec3d& p ){ points.push_back(p); return points.size()-1; };
    inline int addEdge( int a, int b, int type, double l0=-1, double d=0, double mass=0, double kT=0, double kP=0 ){
        TrussEdge edge;
        edge.a=a; edge.b=b; edge.type=type; edge.l0=l0; edge.d=d; edge.mass=mass; edge.kT=kT; edge.kP=kP;
        edges.push_back(edge);
        return edges.size()-1;
    };
    inline int addFace( int a, int b, int c, int type, double area=-1, double mass=0 ){
        TrussFace face;
        face.a=a; face.b=b; face.c=c; face.type=type; face.area=area; face.mass=mass;
        faces.push_back(face);
        return faces.size()-1;
    };

    void clear();
    void sticksFormString( char * str );
    int  loadXYZ( char* fname );
    void affineTransform( Mat3d mat, bool T=false, Vec3d p0=Vec3dZero, Vec3d p=Vec3dZero );

    int pickVertex( const Vec3d &ray0, const Vec3d &hRay ) const;
    int pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R ) const;
    int pickEdge( const Vec3d& ray0, const Vec3d& hRay, double R ) const;

    void panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width );
    void girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width );
    void girder1_caps( int ip0, int ip1, int kind );
    void wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, double width );
    void makeGriders( int nEdges, TrussEdge* edges, Vec3d* points, GirderParams* params, Vec3d * ups );
    void makeGriders( Truss plan, GirderParams* params, Vec3d * ups, std::vector<Vec2i>* ends );
    int makeCylinder( Vec3d p0, Vec3d p1, double r0, double r1, int nPhi=-1, int nL=-1, double dStep=1.0, int edgeTyp=-1, int faceTyp=-1 );
    //void makeCapsular( int n, Vec3d * points, Vec3d * ups, double R, double width, double mass );  // 
    int  addRope( int i, int j, int type, int nsub );

    void updateEdgesLengths();
    void updateFacesAreas();
    void autoBridge(int n, Vec2i * ips, double rmax, int kind );

    Vec2i* getIJs();

    inline Vec2i newBlock(){ Vec2i ps={points.size(),edges.size()}; blocks.push_back( ps ); return ps; };
    void massesToPoints( double* masses );
    // void toMesh( Mesh& mesh );

};

#endif



