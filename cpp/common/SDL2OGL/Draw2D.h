
#ifndef  Draw2D_h
#define  Draw2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

namespace Draw2D{

extern float z_layer; // this can be used for drawing in multiple overlaping layers
                      // should be initialized like this http://stackoverflow.com/questions/19136059/namespace-global-variable-losing-value-c


void drawPoint      ( const Vec2f& vec                   );
void drawPointCross ( const Vec2f& vec, float d          );
void drawVec        ( const Vec2f& vec                   );
void drawVecInPos   ( const Vec2f& v,   const Vec2f& pos );
void drawBody2d     ( const Vec2f& rot, const Vec2f& pos, float l1, float l2 );
void drawLine       ( const Vec2f& p1, const Vec2f& p2   );
void drawTriangle   ( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 );
void drawRectangle  ( float p1x, float p1y, float p2x, float p2y, bool filled );
void drawRectangle  ( const Vec2f& p1, const Vec2f& p2, bool filled );
void drawCircle     ( const Vec2f& center, float radius, int n, bool filled );

void drawPoint_d       ( const Vec2d& vec                   );
void drawPointCross_d  ( const Vec2d& vec, double d         );
void drawVec_d         ( const Vec2d& vec                   );
void drawVecInPos_d    ( const Vec2d& v,   const Vec2d& pos );
void drawBody2d_d      ( const Vec2d& rot, const Vec2d& pos, float l1, float l2 );
void drawLine_d        ( const Vec2d& p1,  const Vec2d& p2  );
void drawTriangle_d    ( const Vec2d& p1,  const Vec2d& p2, const Vec2d& p3 );
void drawRectangle_d   ( const Vec2d& p1,  const Vec2d& p2, bool filled );
void drawCircle_d      ( const Vec2d& center, float radius, int n, bool filled );

void drawPoints        ( int npoints, Vec2d * points );
void drawPoints        ( int npoints, Vec2d * points, float sz );
void drawLines         ( int n, Vec2d * points );
void drawLines         ( int nlinks, int * links, Vec2d * points );
void drawConvexPolygon ( int n, Vec2d * points, bool filled );
void drawPolarFunc     ( double x0, double y0, double fscale, int n, double phi0, double * data );

void plot         ( int n, double * xs, double * ys );
void plot_cross   ( int n, double * xs, double * ys, double sz );

void drawFunc     ( float xmin, float xmax,          int n, Func1d func );
void drawFuncDeriv( float xmin, float xmax, float d, int n, Func1d func );


void drawGrid     ( float xmin, float ymin, float xmax, float ymax, float dx, float dy );

void drawSimplex    ( float x, float y, bool s, float step );
void drawSimplexGrid( int n,  float step );


void drawShape    ( const Vec2d& pos, const Vec2d& rot, int shape );

void renderImage( GLuint itex, const Rect2d& rec );
void drawString ( const char * str, int imin, int imax, float x, float y, float sz, int itex );
void drawString ( const char * str,                     float x, float y, float sz, int itex );


// ==== inline functions

inline void toGLMat( const Vec2d& pos, const Vec2d& rot, float* glMat ){
/*
	glMat[0 ] = (float)rot.ax;   glMat[1 ] = (float)rot.ay;   glMat[2 ] = (float)rot.az;   glMat[3 ]  = 0;
	glMat[4 ] = (float)rot.bx;   glMat[5 ] = (float)rot.by;   glMat[6 ] = (float)rot.bz;   glMat[7 ]  = 0;
	glMat[8 ] = (float)rot.cx;   glMat[9 ] = (float)rot.cy;   glMat[10] = (float)rot.cz;   glMat[11]  = 0;
	glMat[12] = (float)pos. x;   glMat[13] = (float)pos. y;   glMat[14] = (float)pos. z;   glMat[15]  = 1;
*/
	glMat[0 ] = +rot.x; glMat[1 ] = +rot.y; glMat[2 ] = 0;  glMat[3 ] = 0;
	glMat[4 ] = -rot.y; glMat[5 ] = +rot.x; glMat[6 ] = 0;  glMat[7 ] = 0;
	glMat[8 ] = 0;      glMat[9 ] = 0;      glMat[10] = 1;  glMat[11] = 0;
	glMat[12] = pos.x;  glMat[13] = pos.y;  glMat[14] = 0;  glMat[15] = 1;
};

inline void toGLMatCam( const Vec2d& pos, const Vec2d& rot, float* glMat ){
/*
	glMat[0 ] = (float)rot.ax;   glMat[1 ] = (float)rot.bx;   glMat[2 ] = (float)-rot.cx;   glMat[3 ]  = 0;
	glMat[4 ] = (float)rot.ay;   glMat[5 ] = (float)rot.by;   glMat[6 ] = (float)-rot.cy;   glMat[7 ]  = 0;
	glMat[8 ] = (float)rot.az;   glMat[9 ] = (float)rot.bz;   glMat[10] = (float)-rot.cz;   glMat[11]  = 0;
	glMat[12] = (float)-pos. x;  glMat[13] = (float)-pos. y;  glMat[14] = (float)-pos. z;   glMat[15]  = 1;
*/
	glMat[0 ] = +rot.x; glMat[1 ] = -rot.y; glMat[2 ] = 0;  glMat[3 ] = 0;
	glMat[4 ] = +rot.y; glMat[5 ] = +rot.x; glMat[6 ] = 0;  glMat[7 ] = 0;
	glMat[8 ] = 0;      glMat[9 ] = 0;      glMat[10] = 1;  glMat[11] = 0;
	glMat[12] = -pos.x; glMat[13] = -pos.y; glMat[14] = 0;  glMat[15] = 1;
};

}; // namespace Draw2D

#endif

