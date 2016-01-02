
#ifndef  drawMath_h
#define  drawMath_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

// ==== function declarations

void printVec     ( const Vec3d& v );
void printArray   ( int n, double * vec );
void drawPoint    ( const Vec3d& vec );
void drawVec      ( const Vec3d& vec );
void drawVecInPos ( const Vec3d& v, const Vec3d& pos );
void drawLine     ( const Vec3d& p1, const Vec3d& p2 );
void drawVecInPos ( const Vec3f& v, const Vec3f& pos );
void drawTriangle ( const Vec3d& p1, const Vec3d& p2, const Vec3d& p3 );
void drawMatInPos ( const Mat3d& mat, const Vec3d& pos );
//void drawShape    ( const Vec2d& pos, const Vec2d& rot, int shape );
void drawShape    ( const Vec3d& pos, const Mat3d& rot, int shape );

int  drawConeFan        ( int n, float r, const Vec3f& base, const Vec3f& tip );
int  drawCylinderStrip  ( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip );
int  drawSphereTriangle ( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c );
int  drawSphere_oct     ( int n, double r_, const Vec3d& pos_ );	
void drawLines          ( int nlinks, int * links, Vec3d * points );


// ==== inline functions

inline void setColorInt32( uint32_t clr ) {
	const float i255 = 1/255.0f;
	uint8_t b = ( ( clr       ) & 0xFF );   
	uint8_t g = ( ( clr >> 8  ) & 0xFF );
	uint8_t r = ( ( clr >> 16 ) & 0xFF );
	uint8_t a = (   clr >> 24          );
	glColor4f( i255*r, i255*g, i255*b, i255*a );
	//printf( " r %i g %i b %i a %i     %f %f %f %f  \n", r, g, b, a,  i255*r, i255*g, i255*b, i255*a   );
}

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

#endif

