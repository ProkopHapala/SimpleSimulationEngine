
#ifndef  testUtils_h
#define  testUtils_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

void _endl(){ printf("\n"); };

void printArray( int n, double * vec ){	for (int i=0; i<n; i++){ printf( " %f ", vec[i] ); }; printf("\n"); };
void printArray( int n, float  * vec ){	for (int i=0; i<n; i++){ printf( " %f ", vec[i] ); }; printf("\n"); };
void printArray( int n, int    * vec ){	for (int i=0; i<n; i++){ printf( " %i ", vec[i] ); }; printf("\n"); };

void printVec( const Vec2d& v ){ printf( " %f %f \n", v.x, v.y ); };
void printVec( const Vec2f& v ){ printf( " %f %f \n", v.x, v.y );};
void printVec( const Vec2i& v ){ printf( " %i %i \n", v.x, v.y );};

void printVec( const Vec3d& v ){ printf( " %f %f %f \n", v.x, v.y, v.z );};
void printVec( const Vec3f& v ){ printf( " %f %f %f \n", v.x, v.y, v.z );};
void printVec( const Vec3i& v ){ printf( " %i %i %i \n", v.x, v.y, v.z );};

void printQuat( const Quat4d& q ){	printf( " %f %f %f %f \n", q.x, q.y, q.z, q.w ); }
void printQuat( const Quat4f& q ){	printf( " %f %f %f %f \n", q.x, q.y, q.z, q.w ); }
void printQuat( const Quat4i& q ){ 	printf( " %i %i %i %i \n", q.x, q.y, q.z, q.w ); }

void printMat( const Mat3d& mat  ){
	printf( " %f %f %f \n", mat.ax, mat.ay, mat.az );
	printf( " %f %f %f \n", mat.bx, mat.by, mat.bz );
	printf( " %f %f %f \n", mat.cx, mat.cy, mat.cz );
}

_template_T void println(const T& t){t.print(); puts(""); };
_template_T void println(const char* s,const T& t){ printf("%s",s);t.print(); puts(""); };

// CPU ticks timer
// http://stackoverflow.com/questions/6432669/variance-in-rdtsc-overhead

#define DBGL printf( "DEBUG LINE %i \n", __LINE__);

#define DEBUG printf( "DEBUG LINE %i %s %s \n", __LINE__, __FUNCTION__, __FILE__ );

//void dbg(char* s){ printf("DEBUG (%s) \n", s); };

#define DBG(format,args...) { printf("DEBUG "); printf(format,## args); }

/*
int  dbg(int priority, const char *format, ...){
    va_list args;
    va_start(args, format);
    if(priority & PRIO_LOG) vprintf(format, args);
    va_end(args);
}
*/


template<typename Func>
double checkDeriv(Func getEF,const Vec3d p0, double d, Vec3d& fE, Vec3d& f ){
    getEF(p0,f);
    for(int i=0;i<3;i++){
        double E0,E1;
        Vec3d p=p0;
        p .array[i]-=d;   E0=getEF(p,f);
        p .array[i]+=d*2; E1=getEF(p,f);
        p .array[i]-=d;
        fE.array[i]=(E1-E0)/(2*d);
    }
    double err = (f-fE).norm();
    printf( " |f-fE|: %g   fE(%g,%g,%g)   f(%g,%g,%g) \n", err, fE.x, fE.y, fE.z, f.x,f.y,f.z );
    return err;
}


inline uint64_t getCPUticks(){
    uint32_t lo, hi;
    __asm__ __volatile__ (
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "rdtsc\n"
      : "=a" (lo), "=d" (hi)
      :
      : "%ebx", "%ecx" );
    return (uint64_t)hi << 32 | lo;
}

class StopWatch{ public:
    long t1;
    double T;
    void   start(){ t1=getCPUticks(); };
    double stop (){ T=getCPUticks()-t1; return T; };
};

static StopWatch stopWatch;


void genRandomArray( int n, double * vals, double vmin, double vmax ){
    for( int i=0; i<n; i++ ){ vals[i] = randf( vmin, vmax ); }
}

void makeSamples2D( Vec2i ns, Vec3d p0, Vec3d da, Vec3d db, Vec3d *ps ){
    //Vec3d da=a*(1.0d/ns.x);
    //Vec3d db=b*(1.0d/ns.y);
    //printf( "da (%g,%g,%g)\n", da.x,da.y,da.z );
    //printf( "db (%g,%g,%g)\n", db.x,db.y,db.z );
    for(int ib=0; ib<ns.y; ib++){
        Vec3d p = p0+db*ib;
        for(int ia=0; ia<ns.x; ia++){
            *ps = p;
            p.add(da);
            ps++;
        }
    }
}

//template<void func(Vec3d pos, Vec3d& fout)>
template <typename Func>
void sampleFroce(Func func, int n, Vec3d *ps, Vec3d *fs ){
    for(int i=0; i<n; i++){ func(ps[i],fs[i]); }
}

//template<void func(Vec3d pos, Vec3d& fout)>
template <typename Func>
int sampleFroce(Func func, Vec2i ns, Vec3d p0, Vec3d a, Vec3d b, Vec3d*&ps, Vec3d*&fs ){
    int ntot = ns.x*ns.y;
    ps = new Vec3d[ntot];
    fs = new Vec3d[ntot];
    makeSamples2D(ns, p0, a, b, ps);
    //sampleFroce<func>( ntot, ps, fs, func );
    sampleFroce( func, ntot, ps, fs  );
    return ntot;
}

#define SPEED_TEST_FUNC( caption, func, xmin, xmax, ncall ) \
    do{                                    \
    double sum  = 0.0;                     \
    double dx   = (xmax-xmin)/(ncall-1);   \
    long tstart = getCPUticks();           \
    double x = xmin;                       \
    for( int i=0; i<ncall; i++ ){          \
        sum += func ( x );                 \
        x   += dx;                         \
        /*printf("%i %f %f %f \n", i, x, func ( x ), sum  );*/  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \


#define SPEED_TEST_FUNC_ARRAY( caption, func, arr, n, m ) \
    do{                                    \
    double sum  = 0.0;                     \
    long tstart = getCPUticks();           \
    for(int j=0; j<m; j++){                \
        for( int i=0; i<n; i++ ){          \
            sum += func ( arr[i] );        \
            /*printf("%i %f %f %f \n", i, arr[i], func ( arr[i] ), sum  );*/  \
        }                                  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    int ncall = n*m;                       \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \


#define SPEED_TEST_FUNC_NM( caption, func, n, m ) \
    do{                                    \
    double sum  = 0.0;                     \
    long tstart = getCPUticks();           \
    for(int j=0; j<m; j++){                \
        for( int i=0; i<n; i++ ){          \
            sum += func;                   \
        }                                  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    int ncall = n*m;                       \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \



#endif




