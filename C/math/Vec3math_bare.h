
// define these macro befor initialize this
//#define  TYPE  float
//#define  VEC   float3

#include "macros.h"

//inline void VEC_from_array(       VEC* v, const TYPE* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; };
//inline void VEC_to_array  ( const VEC* v, TYPE*       arr ){ arr[0]=v->x; arr[1]=v->y; arr[2]=v->z; };

const VEC METHOD(TPREF,Zero) = (VEC){0,0,0};
const VEC METHOD(TPREF,X)    = (VEC){1,0,0};
const VEC METHOD(TPREF,Y)    = (VEC){0,1,0};
const VEC METHOD(TPREF,Z)    = (VEC){0,0,1};

_inline void METHOD(TPREF, read_array( VEC* v, const TYPE* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; }   )
_inline VEC  METHOD(TPREF, from_array(         const TYPE* arr ){ VEC v; v.x=arr[0]; v.y=arr[1]; v.z=arr[2]; return v; }   )

_inline VEC METHOD(TPREF, clone( VEC  v                ){ VEC o; o.x=v.x; o.y=v.y; o.z=v.z; return o; }   )
_inline VEC METHOD(TPREF, newf ( TYPE f                ){ VEC o; o.x=f; o.y=f; o.z=f; return o; }   )
_inline VEC METHOD(TPREF, new3f( TYPE x, TYPE y, TYPE z){ VEC o; o.x=x; o.y=y; o.z=z; return o; }   )

/*
static inline void METHOD( VEC, setf ( VEC* v, TYPE f                ){ v->x=f; v->y=f; v->z=f; }   )
static inline void METHOD( VEC, set3f( VEC* v, TYPE x, TYPE y, TYPE z){ v->x=x; v->y=y; v->z=z; }   )

static inline void METHOD( VEC, setf( VEC* v, TYPE f ){ v->x=f; v->y=f; v->z=f; }   )
static inline void METHOD( VEC, addf( VEC* v, TYPE f ){ v->x=f; v->y=f; v->z=f; }   )
static inline void METHOD( VEC, mulf( VEC* v, TYPE f ){ v->x=f; v->y=f; v->z=f; }   )

static inline void METHOD( VEC, set3f( VEC* v, TYPE x, TYPE y, TYPE z ){ v->x=x; v->y=y; v->z=z; }   )
static inline void METHOD( VEC, add3f( VEC* v, TYPE x, TYPE y, TYPE z ){ v->x=x; v->y=y; v->z=z; }   )
static inline void METHOD( VEC, mul3f( VEC* v, TYPE x, TYPE y, TYPE z ){ v->x=x; v->y=y; v->z=z; }   )
*/

_inline void METHOD(TPREF, add_p  ( VEC* v, VEC b             ){ v->x+=b.x; v->y+=b.y; v->z+=b.z; }   )
_inline void METHOD(TPREF, sub_p  ( VEC* v, VEC b             ){ v->x-=b.x; v->y-=b.y; v->z-=b.z; }   )
_inline void METHOD(TPREF, mul_p  ( VEC* v, VEC b             ){ v->x*=b.x; v->y*=b.y; v->z*=b.z; }   )

_inline void METHOD(TPREF, addf_p ( VEC* v, TYPE f            ){ v->x+=f; v->y+=f; v->z+=f; }   )
_inline void METHOD(TPREF, mulf_p ( VEC* v, TYPE f            ){ v->x*=f; v->y*=f; v->z*=f; }   )

_inline VEC METHOD(TPREF, add  ( VEC v, VEC b                 ){ v.x+=b.x; v.y+=b.y; v.z+=b.z; return v; }   )
_inline VEC METHOD(TPREF, sub  ( VEC v, VEC b                 ){ v.x-=b.x; v.y-=b.y; v.z-=b.z; return v; }   )
_inline VEC METHOD(TPREF, mul  ( VEC v, VEC b                 ){ v.x*=b.x; v.y*=b.y; v.z*=b.z; return v; }   )

_inline VEC METHOD(TPREF, addf ( VEC v, TYPE f                ){ v.x+=f; v.y+=f; v.z+=f; return v; }   )
_inline VEC METHOD(TPREF, mulf ( VEC v, TYPE f                ){ v.x*=f; v.y*=f; v.z*=f; return v; }   )

_inline void METHOD(TPREF, fma_p ( VEC* v, VEC b, TYPE f      ){ v->x+=b.x*f; v->y+=b.y*f; v->z+=b.z*f; } )

_inline VEC METHOD(TPREF, add3f( VEC v, TYPE x, TYPE y, TYPE z){ v.x+=x; v.y+=y; v.z+=z; return v; }   )
_inline VEC METHOD(TPREF, mul3f( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )

_inline TYPE METHOD(TPREF, norm2( VEC v        ){ return v.x*v.x + v.y*v.y + v.z*v.z; }   )
_inline TYPE METHOD(TPREF, dot  ( VEC a, VEC b ){ return a.x*b.x + a.y*b.y + a.z*b.z; }   )
_inline TYPE METHOD(TPREF, dist2( VEC a, VEC b ){ TYPE dx=a.x-b.x; TYPE dy=a.y-b.y; TYPE dz=a.z-b.z; return dx*dx + dy*dy + dz*dz; }   )

_inline TYPE METHOD(TPREF, norm( VEC v        ){ return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }   )
_inline TYPE METHOD(TPREF, dist( VEC a, VEC b   ){ TYPE dx=a.x-b.x; TYPE dy=a.y-b.y; TYPE dz=a.z-b.z; return sqrt(dx*dx + dy*dy + dz*dz); }   )
_inline TYPE METHOD(TPREF, dirVec( VEC a, VEC b, VEC* o ){ TYPE dx=a.x-b.x; TYPE dy=a.y-b.y; TYPE dz=a.z-b.z; TYPE r=sqrt(dx*dx + dy*dy + dz*dz);   TYPE ir=1/r; o->x=dx*ir; o->y=dy*ir; o->z=dz*ir;   return r; }   )
_inline VEC  METHOD(TPREF, normalized( VEC  v   ){ TYPE r=sqrt(v.x*v.x + v.y*v.y + v.z*v.z);       TYPE ir=1/r; v.x*=ir; v.y*=ir; v.z*=ir;    return v; }   )
_inline TYPE METHOD(TPREF, normalize ( VEC* v   ){ TYPE r=sqrt(v->x*v->x + v->y*v->y + v->z*v->z); TYPE ir=1/r; v->x*=ir; v->y*=ir; v->z*=ir; return r; }   )

_inline VEC METHOD(TPREF, cross( VEC a, VEC b ){ VEC o; o.x=a.y*b.z-a.z*b.y; o.y=a.z*b.x-a.x*b.z; o.z=a.x*b.y-a.y*b.x; return o; }   )


//static inline VEC METHOD(TPREF, norm( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )
//static inline VEC METHOD(TPREF, cos ( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )
//static inline VEC METHOD(TPREF, dist( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )


_inline void METHOD(TPREF, someOrtho(const VEC v, VEC* v1, VEC* v2 ){ )
    TYPE xx = v.x*v.x;
    TYPE yy = v.y*v.y;
    TYPE x,y,z;
    if(xx<yy){
        x =  -yy -v.z*v.z;
        y =  v.x*v.y;
        z =  v.x*v.z;
    }else{
        x =  v.y*v.x;
        y =  -v.z*v.z -xx;
        z =  v.y*v.z;
    }
    v1->x =  x;
    v1->y =  y;
    v1->z =  z;
    v2->x = v.y*z - v.z*y;
    v2->y = v.z*x - v.x*z;
    v2->z = v.x*y - v.y*x;
}

#undef  TYPE
#undef  TPREF
#undef  VEC