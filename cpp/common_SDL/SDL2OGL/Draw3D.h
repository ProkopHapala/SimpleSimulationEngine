
#ifndef  Draw3D_h
#define  Draw3D_h

#include <SDL2/SDL_opengl.h>

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "CMesh.h"
#include "Mesh.h"

#include "Draw.h"

namespace Draw3D{

//typedef Vec3f (*UVfunc)(Vec2f p);

// ==== function declarations

// ========== Float Versions

//void vertex(const Vec3f& v );
//void vertex(const Vec3d& v );
//void color (const Vec3f& v );
//void color (const Vec3d& v );
//void normal(const Vec3f& v );
//void normal(const Vec3d& v );

inline void vertex(const Vec3f& v ){ glVertex3f(v.x,v.y,v.z); }
inline void vertex(const Vec3d& v ){ glVertex3f(v.x,v.y,v.z); }
inline void color (const Vec3f& v ){ glColor3f (v.x,v.y,v.z); }
inline void color (const Vec3d& v ){ glColor3f (v.x,v.y,v.z); }
inline void normal(const Vec3f& v ){ glNormal3f(v.x,v.y,v.z); }
inline void normal(const Vec3d& v ){ glNormal3f(v.x,v.y,v.z); }

void drawPoint     ( const Vec3f& vec                   );
void drawPointCross_bare( const Vec3f& vec, float sz    );
void drawPointCross( const Vec3f& vec, float  sz        );
void drawVec       ( const Vec3f& vec                   );
void drawVecInPos  ( const Vec3f& v,   const Vec3f& pos );
void drawLine      ( const Vec3f& p1,  const Vec3f& p2  );
void drawArrow     ( const Vec3f& p1,  const Vec3f& p2, float sz );

void drawScale     ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& a, float tick, float sza, float szb );

void drawTriangle ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3 );
void drawTriangle ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3, bool filled );

void drawQuad     ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3, const Vec3f& p4, bool filled );


void drawMatInPos ( const Mat3f& mat, const Vec3f& pos, const Vec3f& sc=Vec3fOne );

//void drawShape    ( const Vec3f& pos, const Mat3f&  rot,  int shape, bool transposed = false );
//void drawShape    ( const Vec3f& pos, const Quat4f& qrot, int shape );
//void drawShape    ( const Vec3f& pos, const Quat4f& qrot, const Vec3f& scale, int shape );

void drawShape    ( int shape, const Vec3f& pos, const Mat3f&  rot=Mat3fIdentity, bool transposed = false );
//void drawShape    ( int shape, const Vec3f& pos, const Quat4f& qrot  );
void drawShape    ( int shape, const Vec3f& pos, const Quat4f& qrot, const Vec3f& scale=Vec3fOne );

int  drawCylinderStrip     ( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip );
int  drawCylinderStrip_wire( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip );
int  drawSphereTriangle    ( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c );
int  drawSphereTriangle_wire( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c );

int  drawCircleAxis     ( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R, float dca, float dsa );
int  drawCircleAxis     ( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R );

int  drawConeFan        ( int n, float r,                const Vec3f& base, const Vec3f& tip );
int  drawCone           ( int n, float phi1, float phi2, float r1, float r2, const Vec3f& base, const Vec3f& tip, bool smooth );

//int  drawSphereOctLines ( int n, float R, const Vec3f& pos );
int  drawSphereOctLines ( int n, float R, const Vec3f& pos, const Mat3f& rot=Mat3fIdentity, bool bRGB=false );
int  drawSphere_oct     ( int n, float R, const Vec3f& pos, bool wire=false );
int  drawCapsula        ( Vec3f p0, Vec3f p1,  float r1, float r2, float theta1, float theta2, float dTheta, int nPhi, bool capped );

void drawKite           ( const Vec3f& pos, const Mat3f& rot, float sz );
void drawPanel          ( const Vec3f& pos, const Mat3f& rot, const Vec2f& sz );

void drawText  ( const char * str, const Vec3f& pos, int fontTex, float textSize, int iend=0 );
void drawText3D( const char * str, const Vec3f& pos, const Vec3f& fw, const Vec3f& up, int fontTex, float textSize, int iend=0 );
void drawInt   ( const Vec3d& pos, int i   , int fontTex, float sz=0.02, const char* format="%i\0" );
void drawDouble( const Vec3d& pos, double f, int fontTex, float sz=0.02, const char* format="%g\0" );


void drawBox( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b );
void drawBBox        ( const Vec3f& p0, const Vec3f& p1 );
void drawBBox        ( const Vec3f& p, float r );
void drawTriclinicBox( const Mat3f& lvec_, const Vec3f& c0_, const Vec3f& c1_ );
void drawTriclinicBoxT( const Mat3f& lvec_, const Vec3f& c0_, const Vec3f& c1_ );

void drawAxis( float sc );

// ========== Double Versions

inline void drawPoint     ( const Vec3d& vec                   ){drawPoint((Vec3f)vec); }
inline void drawVec       ( const Vec3d& vec                   ){drawVec  ((Vec3f)vec); }
inline void drawPointCross( const Vec3d& vec, double sz        ){drawPointCross((Vec3f)vec,sz); }
inline void drawVecInPos  ( const Vec3d& v,   const Vec3d& pos ){drawVecInPos((Vec3f)v,(Vec3f)pos); }
inline void drawLine      ( const Vec3d& p1,  const Vec3d& p2  ){drawLine ((Vec3f)p1,(Vec3f)p2); }
inline void drawArrow     ( const Vec3d& p1,  const Vec3d& p2, float sz  ){drawArrow((Vec3f)p1,(Vec3f)p2, sz); }

inline void drawScale     ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& a, double tick, double sza, double szb ){  drawScale( (Vec3f)p1, (Vec3f)p2, (Vec3f)a, tick,sza,szb); };

void vecsInPoss( int n, const Vec3d* vs, const Vec3d* ps, float sc );

inline void drawTriangle ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& p3 ){ drawTriangle( (Vec3f)p1, (Vec3f)p2, (Vec3f)p3 ); };
inline void drawTriangle ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& p3, bool filled ){ drawTriangle( (Vec3f)p1, (Vec3f)p2, (Vec3f)p3, filled ); };
inline void drawTriangle ( const Triangle3D& tri, bool filled ){ drawTriangle( (Vec3f)tri.a, (Vec3f)tri.b, (Vec3f)tri.c, filled ); };
inline void drawQuad     ( const Vec3d& p1,  const Vec3d& p2, const Vec3d& p3, const Vec3d& p4, bool filled ){ drawQuad((Vec3f)p1,(Vec3f)p2,(Vec3f)p3,(Vec3f)p4, filled ); }
inline void drawQuad     ( const Quad3d& qd, bool filled ){ drawQuad(qd.l1.a,qd.l1.b,qd.l2.b,qd.l2.a, filled); }

inline void drawMatInPos ( const Mat3d& mat, const Vec3d& pos, const Vec3d& sc=Vec3dOne ){  drawMatInPos( (Mat3f)mat, (Vec3f)pos, (Vec3f)sc ); };

//inline void drawShape    ( const Vec3d& pos, const Mat3d&  rot,  int shape, bool transposed = false ){ drawShape( (Vec3f)pos, (Mat3f)rot, shape, transposed ); };
//inline void drawShape    ( const Vec3d& pos, const Quat4d& qrot, int shape ){ drawShape( (Vec3f)pos, (Quat4f)qrot, shape); };
//inline void drawShape    ( const Vec3d& pos, const Quat4d& qrot, Vec3d& scale, int shape ){ drawShape( (Vec3f)pos, (Quat4f)qrot, (Vec3f)scale, shape); };

inline void drawShape    ( int shape, const Vec3d& pos, const Mat3d& rot=Mat3dIdentity,  bool transposed = false ){ drawShape(  shape, (Vec3f)pos, (Mat3f)rot,transposed ); };
//inline void drawShape    ( int shape, const Vec3d& pos, const Quat4d& qrot, ){ drawShape( (Vec3f)pos, (Quat4f)qrot, shape); };
inline void drawShape    (  int shape, const Vec3d& pos, const Quat4d& qrot, const Vec3d& scale=Vec3dOne ){ drawShape( shape, (Vec3f)pos, (Quat4f)qrot, (Vec3f)scale); };

void shapeInPoss(  int shape, int n, const Vec3d* pos, const double* sizes=0, const Mat3d& rot=Mat3dIdentity, bool transposed = false );

inline int  drawConeFan        ( int n, float r,                const Vec3d& base,  const Vec3d& tip                                 ){ return drawConeFan( n,             r,      (Vec3f)base, (Vec3f)tip         ); };
inline int  drawCone           ( int n, float phi1, float phi2, float r1, float r2, const Vec3d& base, const Vec3d& tip, bool smooth ){ return drawCone   ( n, phi1, phi2, r1, r2, (Vec3f)base, (Vec3f)tip, smooth ); };

inline int  drawCircleAxis     ( int n, const Vec3d& pos, const Vec3d& v0, const Vec3d& uaxis, double R ){ return drawCircleAxis( n, (Vec3f)pos, (Vec3f)v0, (Vec3f)uaxis, R ); };
//inline int  drawSphereOctLines ( int n, double R, const Vec3d& pos ){ return drawSphereOctLines ( n, R, (Vec3f)pos ); };
inline int  drawSphereOctLines ( int n, double R, const Vec3d& pos, const Mat3d& rot=Mat3dIdentity, bool bRGB=false ){ return drawSphereOctLines ( n, R, (Vec3f)pos, (Mat3f)rot, bRGB ); };
inline int  drawSphere_oct     ( int n, double R, const Vec3d& pos ){ return drawSphere_oct( n, R, (Vec3f)pos ); };

inline void drawKite     ( const Vec3d& pos, const Mat3d& rot, double sz       ){ drawKite ( (Vec3f)pos, (Mat3f)rot, sz ); };
inline void drawPanel    ( const Vec3d& pos, const Mat3d& rot, const Vec2d& sz ){ drawPanel( (Vec3f)pos, (Mat3f)rot, (Vec2f)sz ); };

inline void drawText     ( const char * str, const Vec3d& pos, int fontTex, float textSize, int iend ){ drawText(str, (Vec3f)pos, fontTex, textSize,iend); };

inline void drawBBox         ( const Vec3d& p0, const Vec3d& p1 )                      { drawBBox       ( (Vec3f)p0, (Vec3f)p1 ); };
inline void drawTriclinicBox ( const Mat3d& lvec_, const Vec3d& c0_, const Vec3d& c1_ ){ drawTriclinicBox( (Mat3f)lvec_, (Vec3f)c0_, (Vec3f) c1_ ); };
inline void drawTriclinicBoxT( const Mat3d& lvec_, const Vec3d& c0_, const Vec3d& c1_ ){ drawTriclinicBoxT( (Mat3f)lvec_, (Vec3f)c0_, (Vec3f) c1_ ); };

// ========== Arrays // Not easy to convert

void drawPlanarPolygon( int n, const int * inds, const Vec3d * points );
void drawPolygonNormal( int n, const int * inds, const Vec3d * points );
void drawPolygonBorder( int n, const int * inds, const Vec3d * points );
void drawPlanarPolygon( int ipl, Mesh& mesh );
void drawPolygonBorder( int ipl, Mesh& mesh );
void drawPolygonNormal( int ipl, Mesh& mesh );

//void drawShapeT   ( const Vec3f& pos, const Mat3f&  rot,  int shape );
//void drawShapeT   ( const Vec3d& pos, const Mat3d&  rot,  int shape );
//void drawShapeT   ( const Vec3d& pos, const Quat4d& qrot, int shape );

//int  drawParaboloid     ( Vec3f p0, Vec3f dir, float r, float l, float nR, int nPhi, bool capped );

void drawPolyLine( int n, Vec3d * ps, bool closed=false );

void drawPoints         ( int n, const Vec3d * points, float sz );
void drawLines          ( int nlinks, const int * links, const Vec3d * points );
void drawTriangles      ( int nlinks, const int * links, const Vec3d * point, int mode=0 );
void drawPolygons       ( int nlinks, const int * ns,    const int * links, const Vec3d * points );
void drawMeshWireframe  ( const CMesh& msh );
void drawMeshPolygons   ( const CMesh& msh );
void drawMesh           ( const CMesh& msh, uint32_t cpoly, uint32_t cwire );

void drawVectorArray(int n,const Vec3d* ps,const Vec3d* vs, double sc, double lmax );
void drawScalarArray(int n,const Vec3d* ps,const double* vs, double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );
void drawScalarField(Vec2i ns,const Vec3d* ps,const double* data, double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );

void drawScalarGrid(Vec2i ns,const Vec3d& p0, const Vec3d& a, const Vec3d& b, const double* data,  double vmin=0.0, double vmax=1.0, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );
void drawColorScale( int n,const Vec3d& p0, const Vec3d& Vec3dY, const Vec3d& up=Vec3dX, const uint32_t * colors=&Draw::colors_rainbow[0], int ncol=Draw::ncolors );

void drawSimplexGrid( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs, const double * clrs, int ncolors, const uint32_t * colorscale );
void drawSimplexGridLines( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs );
void drawSimplexGridLinesToned( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs );
void drawRectGridLines( Vec2i n, const Vec3d& p0, const Vec3d& da, const Vec3d& db );

int drawMesh( const Mesh& mesh  );

void drawTetraIso( Vec3f** ps, Quat4d vals );
void drawSimplexLines( Vec3f** ps );

//void drawText( const char * str, const Vec3d& pos, int fontTex, float textSize, int istart, int iend );

//void drawText3D( const char * str, const Vec3d& pos, int fontTex, float textSize, int iend );

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
    //rot.print();
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
        //rot.div(sc);
        rot.mul(sc);
        toGLMat(pos,rot, glMat);
    };
    glMultMatrixf( glMat );
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
inline void drawShapes( int shape, int n, bool transposed, Func func ){
    for(int i=0; i<n; i++){
        Mat3d mat; Vec3d pos;
        func(pos,mat);
        drawShape( shape, pos, mat, transposed );
    }
}

template<typename Func>
inline void drawPBC( const Vec3i& npbc, const Mat3d& lvec, Func func ){
    for(int iz=-npbc.z; iz<=npbc.z; iz++){
        for(int iy=-npbc.y; iy<=npbc.y; iy++){
            for(int ix=-npbc.x; ix<=npbc.x; ix++){
                //builder.pbcShift(ix,iy,iz);
                Vec3d shift = lvec.lincomb( ix, iy, iz );
                glPushMatrix();
                glTranslatef(  shift.x,  shift.y,  shift.z );
                func( );
                //glTranslatef( -shift.x, -shift.y, -shift.z );
                glPopMatrix();
            }
        }
    }
}

}; // namespace Draw3D

#endif

