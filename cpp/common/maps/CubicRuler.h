
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
        pos0  .set( pos0_ );
        mat   .set( mat_ );
        mat.invert_to( invMat );
        printf( "(%3.3f,%3.3f,%3.3f) \n", pos0.x, pos0.y, pos0.z );
        //printf( "(%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", mat.ax,mat.ay,mat.az, mat.bx,mat.by,mat.bz,  mat.cx,mat.cy,mat.cz );
        //printf( "(%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", invMat.ax,invMat.ay,invMat.az, invMat.bx,invMat.by,invMat.bz,  invMat.cx,invMat.cy,invMat.cz );
    }


    Vec3d  hRay,hRayInv,mdRay;
    Vec3i  iRay,diRay;

    void ray_start( const  Vec3d& dirHat, const Vec3d& startPos ){

        Vec3d abc;
        invMat.dot_to( startPos - pos0, abc );

        iRay.a = (int) abc.a;
        iRay.b = (int) abc.b;
        iRay.c = (int) abc.c;

        invMat.dot_to( dirHat,  hRay );
        if( hRay.a < 0 ){ mdRay.a = abc.a-iRay.a; diRay.a=-1; hRay.a=-hRay.a; }else{ mdRay.a = 1 - abc.a + iRay.a; diRay.a=1; };
        if( hRay.b < 0 ){ mdRay.b = abc.b-iRay.b; diRay.b=-1; hRay.b=-hRay.b; }else{ mdRay.b = 1 - abc.b + iRay.b; diRay.b=1; };
        if( hRay.c < 0 ){ mdRay.c = abc.c-iRay.c; diRay.c=-1; hRay.c=-hRay.c; }else{ mdRay.c = 1 - abc.c + iRay.c; diRay.c=1; };

        hRayInv.set_inv( hRay );

    }

   double ray_next( ){
        Vec3d tm;
        tm.set_mul( mdRay, hRayInv );

        if( tm.a < tm.b ){
           if( tm.a < tm.c ){  // a min
                mdRay.a  = 1;
                mdRay.b -= hRay.b*tm.a;
                mdRay.c -= hRay.c*tm.a;
                iRay.a  += diRay.a;
                return tm.a;
           }
        }else{
           if( tm.b < tm.c ){  // b min
                mdRay.a -= hRay.a*tm.b;
                mdRay.b  = 1;
                mdRay.c -= hRay.c*tm.b;
                iRay.b  += diRay.b;
                return tm.b;
           }
        }
        // c min
        mdRay.a -= hRay.a*tm.c;
        mdRay.b -= hRay.b*tm.c;
        mdRay.c = 1;
        iRay.c += diRay.c;
        //hits[i].set_mul( dirHat, t );
        //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], tma, tmb, tmc,       t, hits[i].x, hits[i].y );
        //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], mda, mdb, mdc,       t, hits[i].x, hits[i].y );
        return tm.c;
    }

/*
    double ray_next( ){
        Vec3d tm;
        tm.set_mul( mdRay, hRayInv );

        double dt;

        if( tm.a < tm.b ){
           if( tm.a < tm.c ){  // a min
                dt      += tm.a;
                mdRay.a  = 1;
                mdRay.b -= hRay.b*tm.a;
                mdRay.c -= hRay.c*tm.a;
                iRay.a++;
           }else{             // c min
                dt     = tm.c;
                mdRay.a -= hRay.a*tm.c;
                mdRay.b -= hRay.b*tm.c;
                mdRay.c  = 1;
                iRay.c++;
           }
        }else{
           if( tm.b < tm.c ){  // b min
                dt    = tm.b;
                mdRay.a  -= hRay.a*tm.b;
                mdRay.b = 1;
                mdRay.c -= hRay.c*tm.b;
                iRay.b++;
           }else{              // c min
                dt    = tm.c;
                mdRay.a -= hRay.a*tm.c;
                mdRay.b -= hRay.b*tm.c;
                mdRay.c = 1;
                iRay.c++;
           }
        }
        //hits[i].set_mul( dirHat, t );
        //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], tma, tmb, tmc,       t, hits[i].x, hits[i].y );
        //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], mda, mdb, mdc,       t, hits[i].x, hits[i].y );
        return dt;
    }
*/


/*
    int ray_next( Vec2d dirHat, Vec2d pos0, Vec2d pos1, Vec2d * hits, int * boundaries, int * edges ){
        double t0    = dirHat.dot( pos0 );
        double t1    = dirHat.dot( pos1 );
        double tspan = t1-t0;
        double pa,pb,pc, invPa,invPb,invPc;
        int    ia,ib,ic,i;
        printf( " %f %f \n", step, invStep );
        double mda,mdb,mdc, mta,mtb,mtc;
        pa  = dirHat.dot( { 0.0d        ,1.15470053839*invStep} );
        pb  = dirHat.dot( { 1.0d*invStep,0.57735026919*invStep} );
        pc  = dirHat.dot( {-1.0d*invStep,0.57735026919*invStep} );
        mda = pos0.dot  ( { 0.0d        ,1.15470053839*invStep} );
        mdb = pos0.dot  ( { 1.0d*invStep,0.57735026919*invStep} );
        mdc = pos0.dot  ( {-1.0d*invStep,0.57735026919*invStep} );
        if( pa < 0 ){ pa=-pa; mda = 1-mda; };
        if( pb < 0 ){ pa=-pb; mdb = 1-mdb; };
        if( pc < 0 ){ pc=-pc; mdc = 1-mdc; };
        ia=(int)(mda + MAP_OFFSET);   mda = 1-(mda - (ia - MAP_OFFSET) );
        ib=(int)(mdb + MAP_OFFSET);   mdb = 1-(mdb - (ib - MAP_OFFSET) );
        ic=(int)(mdc + MAP_OFFSET);   mdc = 1-(mdc - (ic - MAP_OFFSET) );
        invPa = 1/pa; invPb = 1/pb; invPc = 1/pc;
        //printf( " t_1,2  %f %f   p_a,b,c %f %f %f  \n", t0, t1, pa, pb, pc );
        printf( " pa invPa \n", pa, invPa );
        double t = 0;
        i=0;
        int ia_,ib_;
        simplexIndexBare( pos0.x, pos0.y, ia_, ib_ );
        while( t<tspan ){
            double tma = mda * invPa;
            double tmb = mdb * invPb;
            double tmc = mdc * invPc;
            //t += tma; boundaries[i] = 0;  mda = 1;
            //t += tmb; boundaries[i] = 1;  mdb = 1;
            //t += tmc; boundaries[i] = 2;  mdc = 1;
            int ii = i<<2;
            if( tma < tmb ){
               if( tma < tmc ){  // a min
                    t    += tma;
                    mda   = 1; mdb -= pb*tma; mdc -= pc*tma;
                    boundaries[i] = 0; ia_++;
                    edges[ii  ] = ia_; edges[ii+1] = ib_;
                    edges[ii+2] = ia_; edges[ii+3] = ib_+1;
               }else{            // c min
                    t    += tmc;
                    mda  -= pa*tmc; mdb -= pb*tmc; mdc = 1;
                    boundaries[i] = 2; ib_++;
                    edges[ii  ] = ia_; edges[ii+1] = ib_;
                    edges[ii+2] = ia_+1; edges[ii+3] = ib_;
               }
            }else{
               if( tmb < tmc ){  // b min
                    t    += tmb;
                    mda  -= pa*tmb; mdb = 1; mdc -= pc*tmb;
                    boundaries[i] = 1;
                    edges[ii  ] = ia_;   edges[ii+1] = ib_+1;
                    edges[ii+2] = ia_+1; edges[ii+3] = ib_;
               }else{            // c min
                    t    += tmc;
                    mda  -= pa*tmc; mdb -= pb*tmc; mdc = 1;
                    boundaries[i] = 2; ib_++;
                    edges[ii  ] = ia_; edges[ii+1] = ib_;
                    edges[ii+2] = ia_+1; edges[ii+3] = ib_;
               }
            }
            hits[i].set( pos0 );
            hits[i].add_mul( dirHat, t );
            //hits[i].set_mul( dirHat, t );
            //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], tma, tmb, tmc,       t, hits[i].x, hits[i].y );
            printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], mda, mdb, mdc,       t, hits[i].x, hits[i].y );
            i++;
        }
        return i;
    }
*/

/*
    while( t<tmax ){
        double tma = mda * invPa;
        double tmb = mdb * invPb;
        double tmc = mdc * invPc;
        int ii = i<<2;
        if( tma < tmb ){
           if( tma < tmc ){  // a min
                t   += tma;
                mda  = 1;
                mdb -= pb*tma;
                mdc -= pc*tma;
                ia++;
           }else{            // c min
                t   += tmc;
                mda -= pa*tmc;
                mdb -= pb*tmc;
                mdc  = 1;
                ic++;
           }
        }else{
           if( tmb < tmc ){  // b min
                t   += tmb;
                mda -= pa*tmb;
                mdb  = 1;
                mdc -= pc*tmb;
                ib++;
           }else{            // c min
                t   += tmc;
                mda -= pa*tmc;
                mdb -= pb*tmc;
                mdc  = 1;
                ic++;
           }
        }
        // do something with ia,ib,ic,mda,mdb,mdc
    }
*/


};

#endif

