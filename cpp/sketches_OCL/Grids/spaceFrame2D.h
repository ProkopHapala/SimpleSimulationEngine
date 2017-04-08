#ifndef  spaceFrame2D_h
#define  spaceFrame2D_h

#include "fastmath.h"
#include "Vec3.h"

//============== Globals

//constexpr int nx   = 4;
//constexpr int ny   = 4;
//constexpr int nx   = 64;
//constexpr int ny   = 64;
constexpr int nx     = 128;
constexpr int ny     = 128;
//constexpr int nx   = 256;
//constexpr int ny   = 256;
//constexpr int nx   = 1024+2;
//constexpr int ny   = 1024+2;
constexpr int ntot = nx*ny;
constexpr int nConstrMax = 16;

int    iConstrains[  nConstrMax];
float   constrains[4*nConstrMax];

float pos   [ntot*4];
float vel   [ntot*4];
float force [ntot*4];
float force_[ntot*4];

float k  = -10.0;
float l0 =  1.0;

inline void addHarmonicForce( const Vec3f& p1, const Vec3f& p2, Vec3f& f, float K ){
    Vec3f d; d.set_sub(p2,p1);
    f.add_mul(d,K);
}

inline void addStickForce( const Vec3f& p1, const Vec3f& p2, Vec3f& f ){
    Vec3f d; d.set_sub(p2,p1);
    float l = d.norm();
    f.add_mul( d, k*(l-l0)/l );
    //printf(" p1 (%g,%g,%g) p2 (%g,%g,%g) l %g f (%g,%g,%g)");
}

void genPos0(int nx, int ny, float * data, float step, float rndStep){
    float scale_val = 1.0/256.0;
    for(int iy=0;iy<ny;iy++){
        for(int ix=0;ix<nx;ix++){
            int i4 = (iy*nx + ix)<<2;
            pos[i4  ] = randf(-rndStep,rndStep) + ix*step;
            pos[i4+1] = randf(-rndStep,rndStep) + iy*step;
            pos[i4+2] = randf(-rndStep,rndStep);
            pos[i4+3] = 1.0f;
        }
    }
}

int findNearest( const Vec3f& p, int n, float * pos ){
    float r2min=1.0e+100f;
    float imin=0;
    for(int i=0;i<n;i++){
        int i4  = i<<2;
        Vec3f d;
        d.set_sub( *((Vec3f*)(pos+i4)), p );
        float r2 = d.norm2();
        if(r2<r2min){
            imin=i;
            r2min=r2;
        }
    }
    return imin;
}

void evalForce(int nx, int ny, float * pos, float * force ){
    const int nx_ = nx-1;
    const int ny_ = ny-1;
    const int idy = nx<<2;
    for(int iy=0;iy<ny;iy++){
        for(int ix=0;ix<nx;ix++){
            int i4      = (iy*nx + ix)<<2;
            float* posi = pos + i4;
            Vec3f f; f.set(0.0f);
            Vec3f p; p = *((Vec3f*)posi);
            if( iy>0   ) addStickForce( *((Vec3f*)(posi-idy)), p, f );
            if( iy<ny_ ) addStickForce( *((Vec3f*)(posi+idy)), p, f );
            if( ix>0   ) addStickForce( *((Vec3f*)(posi-4  )), p, f );
            if( ix<nx_ ) addStickForce( *((Vec3f*)(posi+4  )), p, f );
            f.add(0.0f,-0.05f,0.0f);
            *((Vec3f*)(force + i4)) = f;
        }
    }
}

void move_leapfrog(int n, float * pos, float * vel, float * force, float dt, float damp ){
    for(int i=0;i<n;i++){
        int i4  = i<<2;
        Vec3f v = *((Vec3f*)(vel+i4));
        v.mul(damp);
        v.add_mul( *((Vec3f*)(force+i4)), dt );
        *(Vec3f*)(vel+i4) = v;
        ((Vec3f*)(pos+i4))->add_mul(v, dt);
    }
}

double checkDiff( int n, float * ps, float * p0s ){
    double errSum2 = 0;
    for(int i=0; i<n; i++){
        int i4 = i << 2;
        float err2 = sq(ps[i4]-p0s[i4]) + sq(ps[i4+1]-p0s[i4+1]) + sq(ps[i4+2]-p0s[i4+2]); // + sq(ps[i4+3]-p0s[i4+3]);
        errSum2   += err2;
        if ( err2 > 1e-8 ){
            printf( "%i is (%g,%g,) should be (%g,%g,) \n", i, ps[i4], ps[i4+1], ps[i4+2],    p0s[i4], p0s[i4+1], p0s[i4+2] );
            exit(0);
            break;

        }
    }
    return errSum2;
}


void set_array( int n, float * arr, float f ){ for(int i=0; i<n; i++){ arr[i]=f; } }

//============== Functions

#endif

