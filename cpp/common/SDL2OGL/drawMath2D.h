
#ifndef  drawMath2D_h
#define  drawMath2D_h

#include "Vec2.h"

namespace Draw2D{

void drawPoint      ( const Vec2f& vec                   );
void drawPointCross ( const Vec2f& vec, float d          );
void drawVec        ( const Vec2f& vec                   );
void drawVecInPos   ( const Vec2f& v,   const Vec2f& pos );
void drawLine       ( const Vec2f& p1, const Vec2f& p2   );
void drawTriangle   ( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 );
void drawRectangle  ( const Vec2f& p1, const Vec2f& p2   );
void drawCircle     ( const Vec2f& center, float radius, int n );

void drawPoint_d      ( const Vec2d& vec                   );
void drawPointCross_d ( const Vec2d& vec, double d         );
void drawVec_d        ( const Vec2d& vec                   );
void drawVecInPos_d   ( const Vec2d& v,   const Vec2d& pos );
void drawLine_d       ( const Vec2d& p1,  const Vec2d& p2  );
void drawRectangle_d  ( const Vec2f& p1,  const Vec2f& p2 );
void drawCircle_d     ( const Vec2d& center, float radius, int n );

void drawLines      ( int nlinks, int * links, Vec2d * points );
void drawPolarFunc  ( double x0, double y0, double fscale, int n, double phi0, double * data );

void drawShape    ( const Vec2d& pos, const Vec2d& rot, int shape );


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

