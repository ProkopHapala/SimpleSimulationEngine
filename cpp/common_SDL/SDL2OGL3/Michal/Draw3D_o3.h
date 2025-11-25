
#ifndef  Draw3D_h
#define  Draw3D_h

#include <Renderer.h>



#include <math.h>
#include <cstdlib>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Mesh.h"

#include "Draw.h"

namespace Draw3D{

// ==== function declarations

// ========== Float Versions

inline void vertex(const Vec3f& v ){ opengl1renderer.vertex3f(v.x, v.y, v.z); }
inline void vertex(const Vec3d& v ){ opengl1renderer.vertex3d(v.x, v.y, v.z); }
inline void color (const Vec3f& v ){ opengl1renderer.color3f (v.x, v.y, v.z); }
inline void color (const Vec3d& v ){ opengl1renderer.color3f (v.x, v.y, v.z); }
inline void normal(const Vec3f& v ){ opengl1renderer.normal3f(v.x, v.y, v.z); }
inline void normal(const Vec3d& v ){ opengl1renderer.normal3d(v.x, v.y, v.z); }

void drawPoint     ( const Vec3f& vec                   );
void drawPointCross( const Vec3f& vec, float sz, Vec3f color=Vec3fZero );
void drawVec       ( const Vec3f& vec,           Vec3f color           );
void drawVecInPos  ( const Vec3f& v,   const Vec3f& pos, Vec3f color   );
void drawLine      ( const Vec3f& p1,  const Vec3f& p2,  Vec3f color   );
void drawArrow     ( const Vec3f& p1,  const Vec3f& p2, float sz=0.1 );
void drawTriangle ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3 );
void drawTriangle ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3, bool filled );


void drawMatInPos ( const Mat3f& mat, const Vec3f& pos, const Vec3f& sc=Vec3fOne );

void drawShape    ( int shape, const Vec3f& pos, const Mat3f&  rot=Mat3fIdentity, bool transposed = false );
void drawShape    ( int shape, const Vec3f& pos, const Quat4f& qrot, const Vec3f& scale=Vec3fOne );

int  drawCylinderStrip     ( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip );

int  drawCircleAxis     ( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R, float dca, float dsa );
int  drawCircleAxis     ( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R );

int  drawConeFan        ( int n, float r,                const Vec3f& base, const Vec3f& tip );
int  drawCone           ( int n, float phi1, float phi2, float r1, float r2, const Vec3f& base, const Vec3f& tip, bool smooth );

int  drawSphereOctLines ( int n, float R, const Vec3f& pos, const Mat3f& rot=Mat3fIdentity, bool bRGB=false );
void drawSphereOctLinesInstanced( float r, const std::vector<Vec3f>& ps, Vec3f color=COLOR_WHITE );
void drawSphereOctLinesInstanced( float r, const Vec3d* ps, int n, Vec3f color=COLOR_WHITE );
//void drawSphere_oct     ( int n, float R, const Vec3f& pos, bool wire=false );
void drawSphere ( Vec3f pos, float r, Vec3f color=opengl1renderer.color );

void drawKite           ( const Vec3f& pos, const Mat3f& rot, float sz );
void drawPanel          ( const Vec3f& pos, const Mat3f& rot, const Vec2f& sz );

void drawTextBillboard( const char* str, Vec3f pos, float sz, bool ontop=true, int iend=0 );
inline void drawText  ( const char * str, const Vec3f& pos, int fontTex, float textSize, int iend=0 ) { drawTextBillboard(str, pos, textSize); };
inline void drawText3D( const char * str, const Vec3f& pos, const Vec3f& fw, const Vec3f& up, int fontTex, float textSize, int iend=0 ) { drawTextBillboard(str, pos, textSize); };
void drawInt   ( const Vec3d& pos, int i   , int fontTex, float sz=0.02, const char* format="%i\0" );
void drawDouble( const Vec3d& pos, double f, int fontTex, float sz=0.02, const char* format="%g\0" );

void pointLabels( int n, const Vec3d* ps, int fontTex, float sz=7 );
void drawAxis3D( int n, Vec3d p0, Vec3d dp, double v0, double dval, int fontTex, float tickSz=0.5, float textSz=0.015, const char* format="%g" );
void drawAxis3D( Vec3i ns, Vec3d p0, Vec3d ls, Vec3d v0s, Vec3d dvs, int fontTex, float tickSz=0.5, float textSz=0.015, const char* format="%g" );


void drawBox( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b );
void drawBBox        ( const Vec3f& p0, const Vec3f& p1 );
void drawBBox        ( const Vec3f& p, float r );
void drawTriclinicBox( const Mat3f& lvec_, const Vec3f& c0_=Vec3fZero, const Vec3f& c1_=Vec3fOne );
void drawTriclinicBoxT( const Mat3f& lvec_, const Vec3f& c0_=Vec3fZero, const Vec3f& c1_=Vec3fOne );

void drawAxis( float sc );

// ========== Double Versions

inline void drawPoint     ( const Vec3d& vec                   ){drawPoint((Vec3f)vec); }
inline void drawVec       ( const Vec3d& vec, Vec3f color                ){drawVec((Vec3f)vec,color); }
inline void drawPointCross( const Vec3d& vec, double sz, Vec3f color=Vec3fZero){drawPointCross((Vec3f)vec,sz, color); }
inline void drawVecInPos  ( const Vec3d& v,   const Vec3d& pos, Vec3f color ){drawVecInPos((Vec3f)v,(Vec3f)pos, color); }
inline void drawLine      ( const Vec3d& p1,  const Vec3d& p2, Vec3f color=opengl1renderer.color  ){drawLine ((Vec3f)p1,(Vec3f)p2, color); }
inline void drawArrow     ( const Vec3d& p1,  const Vec3d& p2, float sz=0.1  ){drawArrow((Vec3f)p1,(Vec3f)p2, sz); }

void vecsInPoss( int n, const Vec3d* vs, const Vec3d* ps, float sc, Vec3f color=COLOR_RED );

inline void drawTriangle ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& p3 ){ drawTriangle( (Vec3f)p1, (Vec3f)p2, (Vec3f)p3 ); };
inline void drawTriangle ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& p3, bool filled ){ drawTriangle( (Vec3f)p1, (Vec3f)p2, (Vec3f)p3, filled ); };
inline void drawTriangle ( const Triangle3D& tri, bool filled ){ drawTriangle( (Vec3f)tri.a, (Vec3f)tri.b, (Vec3f)tri.c, filled ); };

inline void drawMatInPos ( const Mat3d& mat, const Vec3d& pos, const Vec3d& sc=Vec3dOne ){  drawMatInPos( (Mat3f)mat, (Vec3f)pos, (Vec3f)sc ); };

inline void drawShape    ( int shape, const Vec3d& pos, const Mat3d& rot=Mat3dIdentity,  bool transposed = false ){ drawShape(  shape, (Vec3f)pos, (Mat3f)rot,transposed ); };
inline void drawShape    (  int shape, const Vec3d& pos, const Quat4d& qrot, const Vec3d& scale=Vec3dOne ){ drawShape( shape, (Vec3f)pos, (Quat4f)qrot, (Vec3f)scale); };

inline int  drawConeFan        ( int n, float r,                const Vec3d& base,  const Vec3d& tip                                 ){ return drawConeFan( n,             r,      (Vec3f)base, (Vec3f)tip         ); };
inline int  drawCone           ( int n, float phi1, float phi2, float r1, float r2, const Vec3d& base, const Vec3d& tip, bool smooth ){ return drawCone   ( n, phi1, phi2, r1, r2, (Vec3f)base, (Vec3f)tip, smooth ); };

inline int  drawCircleAxis     ( int n, const Vec3d& pos, const Vec3d& v0, const Vec3d& uaxis, double R ){ return drawCircleAxis( n, (Vec3f)pos, (Vec3f)v0, (Vec3f)uaxis, R ); };
inline int  drawSphereOctLines ( int n, double R, const Vec3d& pos, const Mat3d& rot=Mat3dIdentity, bool bRGB=false ){ return drawSphereOctLines ( n, R, (Vec3f)pos, (Mat3f)rot, bRGB ); };

inline void drawText     ( const char * str, const Vec3d& pos, int fontTex, float textSize, int iend ){ drawText(str, (Vec3f)pos, fontTex, textSize,iend); };

inline void drawBBox         ( const Vec3d& p0, const Vec3d& p1 )                      { drawBBox       ( (Vec3f)p0, (Vec3f)p1 ); };
inline void drawTriclinicBox ( const Mat3d& lvec_, const Vec3d& c0_=Vec3dZero, const Vec3d& c1_=Vec3dOne ){ drawTriclinicBox( (Mat3f)lvec_, (Vec3f)c0_, (Vec3f) c1_ ); };
inline void drawTriclinicBoxT( const Mat3d& lvec_, const Vec3d& c0_=Vec3dZero, const Vec3d& c1_=Vec3dOne ){ drawTriclinicBoxT( (Mat3f)lvec_, (Vec3f)c0_, (Vec3f) c1_ ); };

// ========== Arrays // Not easy to convert

void drawPlanarPolygon( int n, const int * inds, const Vec3d * points );
void drawPlanarPolygon( int ipl, Mesh& mesh );

void drawPolyLine( int n, Vec3d * ps, bool closed=false );

void drawPoints         ( int n, const Vec3d * points, float sz );
void drawLines          ( int nlinks, const int * links, const Vec3d * points );
void drawTriangles      ( int nlinks, const int * links, const Vec3d * point, int mode=0 );

void drawVectorArray(int n,const Vec3d* ps,const Vec3d*  vs, double sc, double lmax );
void drawVectorArray(int n,const Vec3d* ps,const Quat4f* qs, double sc, double lmax );

void drawScalarArray(int n,const Vec3d* ps,const double* vs, double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );
void drawScalarField(Vec2i ns,const Vec3d*  ps,const double* data, double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );
void drawScalarField(Vec2i ns,const Quat4f* ps,const float* data, int pitch, int offset, double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );

void drawScalarGridLines(Vec2i ns, const Vec3d& p0, const Vec3d& a, const Vec3d& b, const Vec3d& up, const double* data, double sc=1.0, Vec2d vclamp=Vec2d{-1e+300,1e+300} );
void drawScalarGrid(Vec2i ns,const Vec3d& p0, const Vec3d& a, const Vec3d& b, const double* data,  double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );
void drawScalarGrid(Vec2i ns,const Vec3d& p0, const Vec3d& a, const Vec3d& b, const float*  data, int pitch, int offset,  double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );
void drawColorScale( int n,const Vec3d& p0, const Vec3d& Vec3dY, const Vec3d& up=Vec3dX, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );

void drawSimplexGrid( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs, const double * clrs, int ncolors, const uint32_t * colorscale );
void drawSimplexGridLines( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs );
void drawSimplexGridLinesToned( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs );
void drawRectGridLines( Vec2i n, const Vec3d& p0, const Vec3d& da, const Vec3d& db );

void drawTetraIso( Vec3f** ps, Quat4d vals );
void drawCurve    ( float tmin, float tmax,   int n, Func1d3 func );
inline void drawColorScale( int n, Vec3d pos, Vec3d dir, Vec3d up, void (_colorFunc_)(float f) );

// ==== inline functions

// TODO : This could be perhaps moved to Camera

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

inline void toGLMatT( const Vec3f& pos, const Mat3f& rot, float* glMat ){
    //printf("pos (%3.3f,%3.3f,%3.3f)\n", pos.x,pos.y,pos.z);
	glMat[0 ] = rot.ax;   glMat[1 ] = rot.bx;   glMat[2 ] = rot.cx;   glMat[3 ]  = 0;
	glMat[4 ] = rot.ay;   glMat[5 ] = rot.by;   glMat[6 ] = rot.cy;   glMat[7 ]  = 0;
	glMat[8 ] = rot.az;   glMat[9 ] = rot.bz;   glMat[10] = rot.cz;   glMat[11]  = 0;
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

inline void toGLMat( const Vec3f& pos, Mat3f rot, const Vec3f& sc, float* glMat ){
    rot.mul(sc);
    toGLMat( pos, rot, glMat );
};

inline void toGLMatCam( const Vec3f& pos, Mat3f rot, const Vec3f& sc, float* glMat ){
    rot.divT(sc);
    toGLMatCam( pos, rot, glMat );
};

inline void toGLMat   ( const Vec3f& pos, const Quat4f& qrot, float* glMat ){ toGLMat( pos, qrot.toMat(), glMat );    };
inline void toGLMatCam( const Vec3f& pos, const Quat4f& qrot, float* glMat ){ toGLMatCam( pos, qrot.toMat(), glMat ); };

inline void toGLMat   ( const Vec3f& pos, const Quat4f& qrot, const Vec3f& sc, float* glMat ){ toGLMat( pos, qrot.toMat(), sc, glMat );    };
inline void toGLMatCam( const Vec3f& pos, const Quat4f& qrot, const Vec3f& sc, float* glMat ){ toGLMatCam( pos, qrot.toMat(), sc, glMat ); };

inline void rigidTransform( const Vec3f& pos, Mat3f rot, const Vec3f& sc, bool trasposed = false ){
    float glMat[16];
    if( trasposed ){
        rot.mul(sc);
        toGLMatT(pos,rot,glMat);
    }else{
        rot.mul(sc);
        toGLMat(pos,rot, glMat);
    };
    opengl1renderer.multMatrixf( glMat );
}
inline void rigidTransform( const Vec3f& pos, const Quat4f& qrot, const Vec3f& sc, bool trasposed = false ){ rigidTransform( pos, qrot.toMat(), sc, trasposed ); };

inline void toGLMat       ( const Vec3d& pos, const Mat3d& rot,   float* glMat   ){ toGLMat   ( (Vec3f)pos, (Mat3f)rot, glMat ); };
inline void toGLMatT      ( const Vec3d& pos, const Mat3d& rot,   float* glMat   ){ toGLMatT  ( (Vec3f)pos, (Mat3f)rot, glMat ); };
inline void toGLMatCam    ( const Vec3d& pos, const Mat3d& rot,   float* glMat   ){ toGLMatCam( (Vec3f)pos, (Mat3f)rot, glMat ); };
inline void toGLMat       ( const Vec3d& pos, const Quat4d& qrot, float* glMat   ){ toGLMat   ( (Vec3f)pos, ((Quat4f)qrot).toMat(), glMat ); };
inline void toGLMatCam    ( const Vec3d& pos, const Quat4d& qrot, float* glMat   ){ toGLMatCam( (Vec3f)pos, ((Quat4f)qrot).toMat(), glMat ); };
inline void toGLMat       ( const Vec3d& pos, const Mat3d& rot,   const Vec3d& sc, float* glMat   ){ toGLMat   ( (Vec3f)pos, (Mat3f)rot, (Vec3f)sc, glMat ); };
inline void toGLMatCam    ( const Vec3d& pos, const Mat3d& rot,   const Vec3d& sc, float* glMat   ){ toGLMatCam( (Vec3f)pos, (Mat3f)rot, (Vec3f)sc, glMat ); };
inline void toGLMat       ( const Vec3d& pos, const Quat4d& qrot, const Vec3d& sc, float* glMat   ){ toGLMat   ( (Vec3f)pos, ((Quat4f)qrot).toMat(), (Vec3f)sc, glMat ); };
inline void toGLMatCam    ( const Vec3d& pos, const Quat4d& qrot, const Vec3d& sc, float* glMat   ){ toGLMatCam( (Vec3f)pos, ((Quat4f)qrot).toMat(), (Vec3f)sc, glMat ); };
inline void rigidTransform( const Vec3d& pos, const Mat3d& rot,   const Vec3d& sc, bool trasposed = false ){ rigidTransform( (Vec3f)pos, (Mat3f)rot, (Vec3f)sc, trasposed ); };
inline void rigidTransform( const Vec3d& pos, const Quat4d& qrot, const Vec3d& sc, bool trasposed = false ){ rigidTransform( (Vec3f)pos, ((Quat4f)qrot).toMat(), (Vec3f)sc, trasposed ); };

template<typename Func>
inline void drawPBC( const Vec3i& npbc, const Mat3d& lvec, Func func ){
    for(int iz=-npbc.z; iz<=npbc.z; iz++){
        for(int iy=-npbc.y; iy<=npbc.y; iy++){
            for(int ix=-npbc.x; ix<=npbc.x; ix++){
                //builder.pbcShift(ix,iy,iz);
                Vec3d shift = lvec.lincomb( ix, iy, iz );
                opengl1renderer.pushMatrix();
                opengl1renderer.translatef(  shift.x,  shift.y,  shift.z );
                func( {ix,iy,iz} );
                //opengl1renderer.translatef( -shift.x, -shift.y, -shift.z );
                opengl1renderer.popMatrix();
            }
        }
    }
}

template<typename Func>
inline void drawShifts( int n, Vec3d* shifts, int i0, Func func ){
    for(int i=0; i<n; i++){
        const Vec3d shift = shifts[i];
        func( shift, i==i0 );
    }
}

}; // namespace Draw3D

#endif

