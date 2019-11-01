#ifndef  Vec4math_h
#define  Vec4math_h

#include "VectorTypes.h"

#include <math.h>

//static inline void float3_read_array( float3* v, const float* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; };
//static inline float3 float3_from_array( const float* arr ){ float3 v; v.x=arr[0]; v.y=arr[1]; v.z=arr[2]; return v; };

_inline double d_project_beam_to_sphere( double r, double x, double y ){
    double z;
    double r2 = r * r;
    double d2 = x*x + y*y;
    if ( d2 < ( 0.5 * r2 ) ) {
        z = sqrt( r2 - d2 );
    } else {
        double t2 = 0.5 * r;
        z = sqrt( t2 / d2 );
    }
    return z;
}

#define  TYPE float  
#define  TPREF f4  
#define  TPRE3 f3  
#define  VEC  float4  
#define  VEC3 float3
#define  MAT3 float3x3
#include "Vec4math_bare.h"

#define  TYPE double 
#define  TPREF d4
#define  TPRE3 d3
#define  VEC  double4
#define  VEC3 double3 
#define  MAT3 double3x3 
#include "Vec4math_bare.h"

//double3 d3_from_f3( float3  a ){ double3 v; v.x=a.x; v.y=a.y; v.z=a.z; return v; }
//float3  f3_from_d3( double3 a ){ float3  v; v.x=a.x; v.y=a.y; v.z=a.z; return v; }

#endif