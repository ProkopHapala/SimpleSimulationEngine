
#ifndef arrayConvert_h
#define arrayConvert_h

//#include "Vec2.h"
#include "Vec3.h"
//#include ".h"

void float4ToVec3d(int n, float* fs, Vec3d* vs ){
    for(int i=0; i<n; i++ ){  
        int i4 = i<<2;
        vs[i].x=fs[i4+0];
        vs[i].y=fs[i4+1];
        vs[i].y=fs[i4+2];
    }
}

void Vec3dTofloat4(int n, Vec3d* vs, float* fs ){
    for(int i=0; i<n; i++ ){  
        int i4 = i<<2;
        fs[i4+0]=vs[i].x;
        fs[i4+1]=vs[i].y;
        fs[i4+2]=vs[i].z;
    }
}

void Vec3dTofloat8(int n, Vec3d* v1s, Vec3d* v2s, float* fs ){
    for(int i=0; i<n; i++ ){  
        int i8 = i<<3;
        fs[i8+0]=v1s[i].x;
        fs[i8+1]=v1s[i].y;
        fs[i8+2]=v1s[i].z;
        //buff[i8+3]=v1s[i].x;
        fs[i8+4]=v2s[i].x;
        fs[i8+5]=v2s[i].y;
        fs[i8+6]=v2s[i].z;
        //buff[i8+7]=v2s[i].x;
    }
}

#endif

