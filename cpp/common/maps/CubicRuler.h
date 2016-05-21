
#ifndef  CubicRuller_h
#define  CubicRuller_h

#include "Vec3.h"
#include "Mat3.h"





const static int CubicRuler_nFaces  = 6;
const static int CubicRuler_nEdges  = 12;
const static int CubicRuler_nVerts  = 8;

const static int CubicRuler_nNeighs = 1+CubicRuler_nFaces+CubicRuler_nEdges+CubicRuler_nVerts;
constexpr static int  CubicRuler_neighs[ CubicRuler_nNeighs ][3] = {
    {  0,0,0  },
    { -1,0,0  },{ +1,0,0 },  { 0,-1,0 },{ 0,+1,0 },  { 0,0,-1 },{ 0,0,+1 },
    { -1,-1,0 },{ -1,+1,0 }, { +1,-1,0 },{ +1,+1,0 },  { 0,-1,-1 },{ 0,-1,+1 }, { 0,+1,-1 },{ 0,+1,+1 },    { -1,0,-1 },{ +1,0,-1 }, { -1,0,+1 },{ +1,0,+1 },
    { -1,-1,-1  },{ -1,-1,+1 },  { -1,+1,-1 },{ -1,+1,+1 },  { +1,-1,-1  },{ +1,-1,+1 },  { +1,+1,-1 },{ +1,+1,+1 },
};


// ===== point index conversion

inline int_fast64_t xyz2id( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    return ix | ( ((int_fast64_t)iy) << 16 ) | ( ((int_fast64_t)iz) << 32 );
}

inline void id2xyz( int_fast64_t id,  int_fast16_t& ix, int_fast16_t& iy, int_fast16_t& iz ){
    iz=( id & 0xFFFF00000000 ) >> 32;
    iy=( id & 0x0000FFFF0000 ) >> 16;
    ix=( id & 0x00000000FFFF );
}

// ===== class CubicRuler

class CubicRuler {
    public:
    Vec3d pos0;
    //Vec3d step;
    Mat3d mat;
    Mat3d invMat;


    inline void pos2index( const Vec3d& pos, Vec3d& dabc, Vec3i& iabc ) const {
        Vec3d abc;
        invMat.dot_to( pos - pos0, abc );
        iabc.a = (int) abc.a;  dabc.a = abc.a - iabc.a;
        iabc.b = (int) abc.b;  dabc.b = abc.b - iabc.b;
        iabc.c = (int) abc.c;  dabc.c = abc.c - iabc.c;
    }

    inline int_fast64_t pos2index( const Vec3d& pos ) const {
        Vec3d dabc;
        Vec3i iabc;
        pos2index( pos, dabc, iabc );
        //printf( "(%3.3f,%3.3f,%3.3f) (%i,%i,%i)  (%3.3f,%3.3f,%3.3f) %i \n",   pos.x,pos.y,pos.z,  iabc.x,iabc.y,iabc.z,  dabc.x,dabc.y,dabc.z, xyz2id( iabc.x, iabc.y, iabc.z ) );
        return xyz2id( iabc.x, iabc.y, iabc.z );
    }

    inline void index2pos( const Vec3i iabc, const Vec3d& dabc, Vec3d& pos ) const {
        Vec3d abc;
        abc.a = iabc.a + dabc.a;
        abc.b = iabc.b + dabc.b;
        abc.c = iabc.c + dabc.c;
        mat.dot_to( abc, pos );
        pos.add(pos0);
    }

    inline void index2pos( int_fast64_t ind, Vec3d& pos ) const {
        //Vec3i iabc;
        //id2xyz( ind, iabc.x, iabc.y, iabc.z );
        //index2pos( iabc, {0.0,0.0,0.0}, pos );
        int_fast16_t ix,iy,iz;
        id2xyz( ind, ix, iy, iz );
        index2pos( {ix,iy,iz}, {0.0,0.0,0.0}, pos );
    }

    inline void setup( const Vec3d& pos0_, const Mat3d& mat_ ){
        pos0  .set( pos0 );
        mat   .set( mat_ );
        mat.invert_to( invMat );

        //printf( "(%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", mat.ax,mat.ay,mat.az, mat.bx,mat.by,mat.bz,  mat.cx,mat.cy,mat.cz );
        //printf( "(%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", invMat.ax,invMat.ay,invMat.az, invMat.bx,invMat.by,invMat.bz,  invMat.cx,invMat.cy,invMat.cz );
    }

};

#endif

