
#ifndef datatypes2cl_h
#define datatypes2cl_h

#include "datatypes.h"
#include "Vec3.h"
#include "quaternion.h"
#include <CL/cl.h>

inline int add_except(int i,int di,int ino){ return (i==ino)?i:i+di; }
inline Quat4i add_except(Quat4i q, int di, int ino=-1){
    return (Quat4i){ 
        add_except(q.x,di,ino), 
        add_except(q.y,di,ino),
        add_except(q.z,di,ino),
        add_except(q.w,di,ino)
    };
};

void     copy_add(int n, Quat4i* from, Quat4i* to, int i0, int ino=-1 ){ for(int i=0; i<n; i++){ to[i]= add_except( from[i],i0,ino); } }

void     copy    (int n, Quat4i* from, Quat4i* to){ for(int i=0; i<n; i++){ to[i]=from[i]; } }
void     copy    (int n, Quat4f* from, Quat4f* to){ for(int i=0; i<n; i++){ to[i]=from[i]; } }

void     set     (int n, Quat4f* qs, Quat4f v=Quat4fZero ){ for(int i=0; i<n; i++){ qs[i]=v; } }

void     pack    (int n, Quat4d* fs, Quat4f* qs            ){ for(int i=0; i<n; i++){ qs[i]  =(Quat4f)fs[i];            } }
void     pack    (int n, Vec3d*  fs, Quat4f* qs, float K=0 ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];  qs[i].e=K; } }
void   unpack    (int n, Vec3d*  fs, Quat4f* qs            ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f;           } }
double unpack_add(int n, Vec3d*  fs, Quat4f* qs            ){ double E=0; for(int i=0; i<n; i++){ fs[i].add( (Vec3d)qs[i].f ); E+=qs[i].e; }; return E; }

void pack(int n, Vec3d* fs, double* es, Quat4f* qs ){
    for(int i=0; i<n; i++){ 
        float e; if(es){e=es[i];}else{e=0;}
        //float e; if(es){e=es[i];}else{e=0;}
        qs[i].f=(Vec3f)fs[i];
        qs[i].e=e; 
    } 
}
Quat4f* pack(int n, Vec3d* fs, double* es=0){
    Quat4f* qs = new Quat4f[n];
    pack( n, fs, es, qs );
    return qs;
}

void unpack(int n, Vec3d* fs, double* es, Quat4f* qs ){
    for(int i=0; i<n; i++){ 
        const Quat4f& q= qs[i];
        fs[i]      =(Vec3d)q.f;
        if(es)es[i]=       q.e; 
    } 
}

// coversions between OpenCL and C++ types
void v2f4  ( const Vec3d& v, float4& f4 ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };                      // pack Vec3d to float4
void f4toV3( const cl_float4& f4, Vec3d& v ){ v.x=f4.s[0]; v.y=f4.s[1]; v.z=f4.s[2]; };                               // unpack float4 to Vec3d
void v2i4  ( const Vec3i& v, int4& f4   ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };                        // pack Vec3i to int4
//void v2f4( const Vec3d& v, cl_float4& f4 ){ f4.s[0]=(cl_float)v.x; f4.s[1]=(cl_float)v.y; f4.s[2]=(cl_float)v.z; }; // pack Vec3d to float4
cl_float4 cl_f4( const Vec3d& v ){ return (cl_float4){(cl_float)v.x,(cl_float)v.y,(cl_float)v.z,0.f}; };              // pack Vec3d to float4

// see https://stackoverflow.com/questions/33119233/opencl-using-struct-as-kernel-argument
typedef struct __attribute__ ((packed)) cl_Mat3{
    cl_float4 a;
    cl_float4 b;
    cl_float4 c;
}st_foo;

inline void printMat( const cl_Mat3& mat  ){
	printf( " %f %f %f \n", mat.a.s[0], mat.a.s[1], mat.a.s[2] );
	printf( " %f %f %f \n", mat.b.s[0], mat.b.s[1], mat.b.s[2] );
	printf( " %f %f %f \n", mat.c.s[0], mat.c.s[1], mat.c.s[2] );
}

// pack Mat3d to OpenCL cl_Mat3
void Mat3_to_cl( const Mat3d& m, cl_Mat3& clm ){ 
    //v2f4( m.a, clm.a ); v2f4( m.b, clm.b ); v2f4( m.c, clm.c ); 
    clm.a = cl_f4( m.a );
    clm.b = cl_f4( m.b );
    clm.c = cl_f4( m.c );
}

// unpack OpenCL cl_Mat3 to Mat3d
void Mat3_from_cl( Mat3d& m, const cl_Mat3& clm ){ 
    f4toV3( clm.a, m.a );
    f4toV3( clm.b, m.b );
    f4toV3( clm.c, m.c );
}



#endif
