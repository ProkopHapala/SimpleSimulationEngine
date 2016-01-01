
#ifndef  drawMath2D_h
#define  drawMath2D_h

#include "Vec2.h"

void printVec       ( const Vec2d& v );
void drawPoint      ( const Vec2f& vec );
void drawPointCross ( const Vec2f& vec, float d );
void drawVec        ( const Vec2f& vec );
void drawVecInPos   ( const Vec2f& v, const Vec2f& pos );
void drawLine       ( const Vec2f& p1, const Vec2f& p2 );
void drawTriangle   ( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 );

void drawPoint      ( const Vec2d& vec                   );
void drawPointCross ( const Vec2d& vec, double d         );
void drawVec        ( const Vec2d& vec                   );
void drawVecInPos   ( const Vec2d& v,   const Vec2d& pos );
void drawLine       ( const Vec2d& p1,  const Vec2d& p2  );

void drawLines      ( int nlinks, int * links, Vec2d * points );
void drawPolarFunc  ( double x0, double y0, double fscale, int n, double phi0, double * data );

#endif

