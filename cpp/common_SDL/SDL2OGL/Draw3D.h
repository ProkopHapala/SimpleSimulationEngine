
#ifndef  Draw3D_h
#define  Draw3D_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Mesh.h"

namespace Draw3D{

// ==== function declarations

void drawPoint     ( const Vec3d& vec                   );
void drawPointCross( const Vec3f& vec, float  sz        );
void drawPointCross( const Vec3d& vec, double sz        );
void drawVec       ( const Vec3d& vec                   );
void drawVecInPos  ( const Vec3d& v,   const Vec3d& pos );
void drawVecInPos  ( const Vec3f& v,   const Vec3f& pos );
void drawLine      ( const Vec3d& p1,  const Vec3d& p2  );
void drawScale     ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& a, double tick, double sza, double szb );

void drawTriangle ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& p3 );
void drawPlanarPolygon( int n, const int * inds, const Vec3d * points );
void drawMatInPos ( const Mat3d& mat, const Vec3d& pos );
//void drawShape    ( const Vec2d& pos, const Vec2d& rot, int shape );
void drawShape    ( const Vec3f& pos, const Mat3f& rot, int shape );
void drawShape    ( const Vec3d& pos, const Mat3d& rot, int shape );
void drawShape    ( const Vec3d& pos, const Quat4d& qrot, int shape );

int  drawConeFan        ( int n, float r,                const Vec3f& base, const Vec3f& tip );
int  drawCone           ( int n, float phi1, float phi2, float r1, float r2, const Vec3f& base, const Vec3f& tip, bool smooth );

int  drawCylinderStrip  ( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip );
int  drawSphereTriangle ( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c );
int  drawSphere_oct     ( int n, double r_, const Vec3d& pos_ );

int drawCircleAxis      ( int n, const Vec3d& pos, const Vec3d& v0, const Vec3d& uaxis, double dca, double dsa );
int drawCircleAxis      ( int n, const Vec3d& pos, const Vec3d& v0, const Vec3d& uaxis );
int drawSphereOctLines  ( int n, double r, const Vec3d& pos );

void drawPolyLine( int n, Vec3d * ps, bool closed=false );

void drawLines          ( int nlinks, const int * links, const Vec3d * points );
void drawTriangles      ( int nlinks, const int * links, const Vec3d * points );
void drawPolygons       ( int nlinks, const int * ns,    const int * links, const Vec3d * points );

void drawKite           ( const Vec3d& pos, const Mat3d& rot, double sz );

void drawSimplexGrid( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs, const double * clrs, int ncolors, const uint32_t * colorscale );
void drawSimplexGridLines( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs );

void drawRectGridLines( Vec2i n, const Vec3d& p0, const Vec3d& da, const Vec3d& db );

int drawMesh( const Mesh& mesh  );

void drawText( const char * str, const Vec3d& pos, int fontTex, float textSize, int istart, int iend );

// from drawUtils.h
void drawBox( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b );
void drawBBox( const Vec3d& p0, const Vec3d& p1 );
int  makeBoxList( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b  );
void drawAxis( float sc );


void drawCurve    ( float tmin, float tmax,          int n, Func1d3 func );

// ==== inline functions

inline void drawColorScale( int n, Vec3d pos, Vec3d dir, Vec3d up, void (_colorFunc_)(float f) ){
    glBegin(GL_TRIANGLE_STRIP);
    double d = 1.0/(n-1);
    for(int i=0; i<n; i++){
        double f = i*d;
        _colorFunc_( f );
        //glColor3f(1.0,1.0,1.0);
        Vec3d p = pos + dir*f;
        glVertex3f( (float)(p.x     ),(float)( p.y     ),(float)( p.z     ) );
        glVertex3f( (float)(p.x+up.x),(float)( p.y+up.y),(float)( p.z+up.z) );
        //printf( "(%g,%g,%g) (%g,%g,%g) \n", p.x, p.y, p.z, (float)(pos.x+up.x),(float)( pos.y+up.y),(float)( pos.z+up.z)  );
    }
    glEnd();
}


/*
inline void toGLMat( const Vec3d& pos, const Mat3d& rot, float* glMat ){
	glMat[0 ] = (float)rot.ax;   glMat[1 ] = (float)rot.ay;   glMat[2 ] = (float)rot.az;   glMat[3 ]  = 0;
	glMat[4 ] = (float)rot.bx;   glMat[5 ] = (float)rot.by;   glMat[6 ] = (float)rot.bz;   glMat[7 ]  = 0;
	glMat[8 ] = (float)rot.cx;   glMat[9 ] = (float)rot.cy;   glMat[10] = (float)rot.cz;   glMat[11]  = 0;
	glMat[12] = (float)pos. x;   glMat[13] = (float)pos. y;   glMat[14] = (float)pos. z;   glMat[15]  = 1;
};

inline void toGLMatCam( const Vec3d& pos, const Mat3d& rot, float* glMat ){
	glMat[0 ] = (float)rot.ax;   glMat[1 ] = (float)rot.bx;   glMat[2 ] = (float)-rot.cx;   glMat[3 ]  = 0;
	glMat[4 ] = (float)rot.ay;   glMat[5 ] = (float)rot.by;   glMat[6 ] = (float)-rot.cy;   glMat[7 ]  = 0;
	glMat[8 ] = (float)rot.az;   glMat[9 ] = (float)rot.bz;   glMat[10] = (float)-rot.cz;   glMat[11]  = 0;
	glMat[12] = (float)-pos. x;  glMat[13] = (float)-pos. y;  glMat[14] = (float)-pos. z;   glMat[15]  = 1;
};
*/

inline void toGLMat( const Vec3f& pos, const Mat3f& rot, float* glMat ){
    //printf("pos (%3.3f,%3.3f,%3.3f)\n", pos.x,pos.y,pos.z);
	glMat[0 ] = rot.ax;   glMat[1 ] = rot.ay;   glMat[2 ] = rot.az;   glMat[3 ]  = 0;
	glMat[4 ] = rot.bx;   glMat[5 ] = rot.by;   glMat[6 ] = rot.bz;   glMat[7 ]  = 0;
	glMat[8 ] = rot.cx;   glMat[9 ] = rot.cy;   glMat[10] = rot.cz;   glMat[11]  = 0;
	glMat[12] = pos. x;   glMat[13] = pos. y;   glMat[14] = pos. z;   glMat[15]  = 1;
    //glMat[0 ] = rot.ax;   glMat[1 ] = rot.ay;   glMat[2 ] = rot.az;   glMat[3 ]  = pos.x;
	//glMat[4 ] = rot.bx;   glMat[5 ] = rot.by;   glMat[6 ] = rot.bz;   glMat[7 ]  = pos.y;
	//glMat[8 ] = rot.cx;   glMat[9 ] = rot.cy;   glMat[10] = rot.cz;   glMat[11]  = pos.z;
	//glMat[12] = 0;        glMat[13] = 0;        glMat[14] = 0;        glMat[15]  = 1;
};

inline void toGLMatCam( const Vec3f& pos, const Mat3f& rot, float* glMat ){
	glMat[0 ] = rot.ax;   glMat[1 ] = rot.bx;   glMat[2 ] = -rot.cx;   glMat[3 ]  = 0;
	glMat[4 ] = rot.ay;   glMat[5 ] = rot.by;   glMat[6 ] = -rot.cy;   glMat[7 ]  = 0;
	glMat[8 ] = rot.az;   glMat[9 ] = rot.bz;   glMat[10] = -rot.cz;   glMat[11]  = 0;
	glMat[12] = -pos. x;  glMat[13] = -pos. y;  glMat[14] = -pos. z;   glMat[15]  = 1;
};


inline void toGLMat( const Vec3d& pos, const Mat3d& rot, float* glMat ){
    Vec3f pos_; convert( pos, pos_ );
    Mat3f rot_; convert( rot, rot_ );
    toGLMat( pos_, rot_, glMat );
};

inline void toGLMatCam( const Vec3d& pos, const Mat3d& rot, float* glMat ){
    Vec3f pos_; convert( pos, pos_ );
    Mat3f rot_; convert( rot, rot_ );
    toGLMatCam( pos_, rot_, glMat );
};

inline void toGLMat( const Vec3d& pos, const Quat4d& qrot, float* glMat ){
    Vec3f pos_;    convert( pos, pos_ );
    Quat4f qrot_;  convert( qrot, qrot_ );
    Mat3f  rot_;   qrot_.toMatrix_unitary2( rot_ );
    toGLMat( pos_, rot_, glMat );
};

inline void toGLMatCam( const Vec3d& pos, const Quat4d& qrot, float* glMat ){
    Vec3f pos_;    convert( pos, pos_ );
    Quat4f qrot_;  convert( qrot, qrot_ );
    Mat3f  rot_;   qrot_.toMatrix_unitary2( rot_ );
    toGLMatCam( pos_, rot_, glMat );
};

inline void toGLMat( const Vec3d& pos, const Mat3d& rot,  const Vec3d& sc, float* glMat ){
    Vec3f pos_; convert( pos, pos_ );
    Vec3f sc_;  convert( sc, sc_ );
    Mat3f rot_; convert( rot, rot_ );
    rot_.mul(sc_);
    toGLMat( pos_, rot_, glMat );
};

inline void toGLMatCam( const Vec3d& pos, const Mat3d& rot,  const Vec3d& sc, float* glMat ){
    Vec3f pos_; convert( pos, pos_ );
    Vec3f sc_;  convert( sc, sc_ );
    Mat3f rot_; convert( rot, rot_ );
    sc_.set_inv(sc_);
    rot_.mulT(sc_);
    toGLMatCam( pos_, rot_, glMat );
};


}; // namespace Draw3D

#endif

