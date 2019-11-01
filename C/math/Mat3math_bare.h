
// define these macro befor initialize this
//#define  TYPE  float
//#define  VEC   float3

#include "macros.h"

const MAT METHOD(TPREF,Zero) = (MAT){0,0,0, 0,0,0, 0,0,0};
const MAT METHOD(TPREF,Eye ) = (MAT){1,0,0, 0,1,0, 0,0,1};

//const MAT METHOD_NAME(TPREF,X)    = (MAT){1,0,0};
//const MAT METHOD_NAME(TPREF,Y)    = (MAT){0,1,0};
//const MAT METHOD_NAME(TPREF,Z)    = (MAT){0,0,1};

/*

_inline void METHOD(TPREF, read_array( VEC* v, const TYPE* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; }   )
_inline VEC  METHOD(TPREF, from_array(         const TYPE* arr ){ VEC v; v.x=arr[0]; v.y=arr[1]; v.z=arr[2]; return v; }   )

_inline VEC METHOD(TPREF, clone( VEC  v                ){ VEC o; o.x=v.x; o.y=v.y; o.z=v.z; return o; }   )
_inline VEC METHOD(TPREF, newf ( TYPE f                ){ VEC o; o.x=f; o.y=f; o.z=f; return o; }   )
_inline VEC METHOD(TPREF, new3f( TYPE x, TYPE y, TYPE z){ VEC o; o.x=x; o.y=y; o.z=z; return o; }   )


_inline VEC METHOD(TPREF, add  ( VEC v, VEC b                 ){ v.x+=b.x; v.y+=b.y; v.z+=b.z; return v; }   )
_inline VEC METHOD(TPREF, mul  ( VEC v, VEC b                 ){ v.x*=b.x; v.y*=b.y; v.z*=b.z; return v; }   )

_inline VEC METHOD(TPREF, addf ( VEC v, TYPE f                ){ v.x+=f; v.y+=f; v.z+=f; return v; }   )
_inline VEC METHOD(TPREF, mulf ( VEC v, TYPE f                ){ v.x*=f; v.y*=f; v.z*=f; return v; }   )

_inline VEC METHOD(TPREF, add3f( VEC v, TYPE x, TYPE y, TYPE z){ v.x+=x; v.y+=y; v.z+=z; return v; }   )
_inline VEC METHOD(TPREF, mul3f( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )

_inline TYPE METHOD(TPREF, norm2( VEC v        ){ return v.x*v.x + v.y*v.y + v.z*v.z; }   )
_inline TYPE METHOD(TPREF, dot  ( VEC a, VEC b ){ return a.x*b.x + a.y*b.y + a.z*b.z; }   )
_inline TYPE METHOD(TPREF, dist2( VEC a, VEC b ){ TYPE dx=a.x-b.x; TYPE dy=a.y-b.y; TYPE dz=a.z-b.z; return dx*dx + dy*dy + dz*dz; }   )

_inline VEC METHOD(TPREF, cross( VEC a, VEC b ){ VEC o; o.x=a.y*b.z-a.z*b.y; o.y=a.z*b.x-a.x*b.z; o.z=a.x*b.y-a.y*b.x; return o; }   )

*/

#undef  TPREF
#undef  TYPE
#undef  VEC
#undef  MAT