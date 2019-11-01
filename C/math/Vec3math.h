#ifndef  Vec3math_h
#define  Vec3math_h

#include "VectorTypes.h"
#include <math.h>

//static inline void float3_read_array( float3* v, const float* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; };
//static inline float3 float3_from_array( const float* arr ){ float3 v; v.x=arr[0]; v.y=arr[1]; v.z=arr[2]; return v; };

#define  TYPE float  
#define  TPREF f3  
#define  VEC  float3 
#include "Vec3math_bare.h"

#define  TYPE double 
#define  TPREF d3
#define  VEC  double3 
#include "Vec3math_bare.h"

double3 d3_from_f3( float3  a ){ double3 v; v.x=a.x; v.y=a.y; v.z=a.z; return v; }
float3  f3_from_d3( double3 a ){ float3  v; v.x=a.x; v.y=a.y; v.z=a.z; return v; }

#endif