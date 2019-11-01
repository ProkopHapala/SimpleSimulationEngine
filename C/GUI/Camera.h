
//#include "Vec3math.h"
//#include "fastMath.h"

#ifndef  Camera_h
#define  Camera_h

#include "VectorTypes.h"

typedef struct Camera{
    float3    pos; // = (float3){0.0f,0.0f,-50.0f};
    float3x3  rot; //= Mat3fIdentity;
    float  zoom  ; //= 10.0f;
    float  aspect; //= 1.0;
    float  zmin  ; //= 10.0;
    float  zmax  ; //= 10000.0;
} Camera;




#endif