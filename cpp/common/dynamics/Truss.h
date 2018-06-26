
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
    int n=5;
    int m=2;
    int kind=0;
    double widthX = 0.5;
    double widthY = 0.5;
    int kind_long   = 0;
    int kind_perp   = 1;
    int kind_zigIn  = 2;
    int kind_zigOut = 3;
};

class TrussEdge{ public:
    int a,b;
    int type;

    void fromString(char* str){ sscanf ( str, "%i %i %i",&a,&b,&type); a--; b--; }
    void toString  (char* str){ sprintf( str, "%i %i %i", a+1, b+1, type); }
    void print     (         ){  printf(      "%i %i %i", a+1, b+1, type); }
};

class Truss{ public:
    int    rendered_shape;
    std::vector<Vec3d>        points;
    std::vector<TrussEdge>     edges;
    std::vector<Vec2i>        blocks;
    std::vector<int>          removed_points;
    std::vector<int>          removed_edges;

    void clear();
    void sticksFormString( char * str );
    int loadXYZ( char* fname );
    void affineTransform( Mat3d mat, bool T );

    int pickVertex( const Vec3d &ray0, const Vec3d &hRay ) const;
    int pickVertex( const Vec3d& ray0, const Vec3d& hRay, double R ) const;
    int pickEdge( const Vec3d& ray0, const Vec3d& hRay, double R ) const;

    void panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width );
    void girder1( Vec3d p0, Vec3d p1, Vec3d up, int n, double width );
    void girder1_caps( int ip0, int ip1, int kind );
    void wheel( Vec3d p0, Vec3d p1, Vec3d ax, int n, double width );
    void makeGriders( int nEdges, TrussEdge* edges, Vec3d* points, GirderParams* params, Vec3d * ups );
    void makeGriders( Truss plan, GirderParams* params, Vec3d * ups, std::vector<Vec2i>* ends );
    void autoBridge(int n, Vec2i * ips, double rmax, int kind );


    Vec2i* getIJs();


    inline Vec2i newBlock(){ Vec2i ps={points.size(),edges.size()}; blocks.push_back( ps ); return ps; };

};

#endif



