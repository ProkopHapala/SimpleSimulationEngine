#include "macros.h"

#define TRACKBALLSIZE ( 0.8 )


//inline void VEC_from_array(       VEC* v, const TYPE* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; };
//inline void VEC_to_array  ( const VEC* v, TYPE*       arr ){ arr[0]=v->x; arr[1]=v->y; arr[2]=v->z; };

const VEC METHOD(TPREF,Zero) = (VEC){0,0,0,0};
const VEC METHOD(TPREF,X)    = (VEC){1,0,0,0};
const VEC METHOD(TPREF,Y)    = (VEC){0,1,0,0};
const VEC METHOD(TPREF,Z)    = (VEC){0,0,1,0};
const VEC METHOD(TPREF,W)    = (VEC){0,0,0,1};
const VEC METHOD(TPREF,Eye)  = (VEC){0,0,0,1};

_inline void METHOD(TPREF, read_array( VEC* v, const TYPE* arr ){ v->x=arr[0]; v->y=arr[1]; v->z=arr[2]; v->w=arr[3]; }   )
_inline VEC  METHOD(TPREF, from_array(         const TYPE* arr ){ VEC v; v.x=arr[0]; v.y=arr[1]; v.z=arr[2]; v.w=arr[3]; return v; }   )

_inline VEC METHOD(TPREF, clone( VEC  v                        ){ VEC o; o.x=v.x; o.y=v.y; o.z=v.z; o.w=v.w;  return o; }   )
_inline VEC METHOD(TPREF, newf ( TYPE f                        ){ VEC o; o.x=f; o.y=f; o.z=f; o.w=f; return o; }   )
_inline VEC METHOD(TPREF, new3f( TYPE x, TYPE y, TYPE z, TYPE w){ VEC o; o.x=x; o.y=y; o.z=z; o.w=w; return o; }   )


_inline VEC METHOD(TPREF, add  ( VEC v, VEC b                 ){ v.x+=b.x; v.y+=b.y; v.z+=b.z; v.w+=b.w; return v; }   )
_inline VEC METHOD(TPREF, mul  ( VEC v, VEC b                 ){ v.x*=b.x; v.y*=b.y; v.z*=b.z; v.w*=b.w; return v; }   )

_inline VEC METHOD(TPREF, addf ( VEC v, TYPE f                ){ v.x+=f; v.y+=f; v.z+=f; v.w+=f; return v; }   )
_inline VEC METHOD(TPREF, mulf ( VEC v, TYPE f                ){ v.x*=f; v.y*=f; v.z*=f; v.w*=f; return v; }   )

_inline VEC METHOD(TPREF, add3f( VEC v, TYPE x, TYPE y, TYPE z,TYPE w){ v.x+=x; v.y+=y; v.z+=z; return v; }   )
_inline VEC METHOD(TPREF, mul3f( VEC v, TYPE x, TYPE y, TYPE z,TYPE w){ v.x*=x; v.y*=y; v.z*=z; return v; }   )

_inline TYPE METHOD(TPREF, norm2( VEC v        ){ return v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w; }   )
_inline TYPE METHOD(TPREF, dot  ( VEC a, VEC b ){ return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }   )
_inline TYPE METHOD(TPREF, dist2( VEC a, VEC b ){ TYPE dx=a.x-b.x; TYPE dy=a.y-b.y; TYPE dz=a.z-b.z; TYPE dw=a.w-b.w; return dx*dx + dy*dy + dz*dz + dw*dw; }   )



// ================= Quaternion

_inline VEC METHOD(TPREF, fromCosAngleAxis( TYPE scal_prod, const VEC3 axis ){ )
// we assume -phi instead of phi!!!, minus effectively implies sa -> -s
    VEC o;
    TYPE cos_cutoff = 1 - 1e-6;
    //TYPE cosphi, sinphi, sa, phi, sgn_sinphi;
    TYPE cosphi, sgn_sinphi, sa;
    TYPE ir    = 1.0 / sqrt( METHOD(TPRE3, norm2(axis) ) );
    VEC3  hat  = METHOD(TPRE3, mulf( axis , ir ) );
    cosphi     = scal_prod;
    sgn_sinphi = 1.0; // ?
    if( cosphi > cos_cutoff ){
        sa = 0;  o.w = 1;
    } else if( cosphi < -( cos_cutoff ) ){
        sa = -1; o.w = 0;
    } else {
        sa = + sqrt( (1 - cosphi) / 2.0 );
        o.w  = - sqrt( (1 + cosphi) / 2.0 ) * sgn_sinphi;
//			sa = -sa; w = -w;
    }
    //VEC o;
    o.x = sa * hat.x;
    o.y = sa * hat.y;
    o.z = sa * hat.z;
    return o;
}


_inline MAT3 METHOD(TPREF, toMat3( VEC v ){ )
    TYPE r2 = v.w*v.w + v.x*v.x + v.y*v.y + v.z*v.z;
    //T s  = (r2 > 0) ? 2d / r2 : 0;
    TYPE s  = 2 / r2;
    // compute xs/ys/zs first to save 6 multiplications, since xs/ys/zs
    // will be used 2-4 times each.
    TYPE xs = v.x * s;  TYPE ys = v.y * s;  TYPE zs = v.z * s;
    TYPE xx = v.x * xs; TYPE xy = v.x * ys; TYPE xz = v.x * zs;
    TYPE xw = v.w * xs; TYPE yy = v.y * ys; TYPE yz = v.y * zs;
    TYPE yw = v.w * ys; TYPE zz = v.z * zs; TYPE zw = v.w * zs;
    // using s=2/norm (instead of 1/norm) saves 9 multiplications by 2 here
    MAT3 m;
    m.xx = 1 - (yy + zz);
    m.xy =     (xy - zw);
    m.xz =     (xz + yw);
    m.yx =     (xy + zw);
    m.yy = 1 - (xx + zz);
    m.yz =     (yz - xw);
    m.zx =     (xz - yw);
    m.zy =     (yz + xw);
    m.zz = 1 - (xx + yy);
    return m;
}

_inline VEC METHOD(TPREF, qmul( const VEC a, const VEC b ){ )
    VEC o;
    o.x =  a.x * b.w + a.y * b.z - a.z * b.y + a.w * b.x;
    o.y = -a.x * b.z + a.y * b.w + a.z * b.x + a.w * b.y;
    o.z =  a.x * b.y - a.y * b.x + a.z * b.w + a.w * b.z;
    o.w = -a.x * b.x - a.y * b.y - a.z * b.z + a.w * b.w;
    return o;
}

_inline VEC METHOD(TPREF, fromTrackball( TYPE p1x, TYPE p1y, TYPE p2x, TYPE p2y ){ )
    //VEC3  axis; // axis of rotation
    //T phi;  // angle of rotation
    //VEC3  p1, p2, d;
    //T t;
    //if( p1x == p2x && p1y == p2y ){	setOne(); return; }
    TYPE dx=p2x-p1x;
    TYPE dy=p2y-p1y;
    if( ( dx*dx + dy*dy ) < 1e-8 ){ return METHOD(TPREF, Eye); }
    VEC3 p1   = (VEC3){ p1x, p1y, d_project_beam_to_sphere( TRACKBALLSIZE, p1x, p1y ) };
    VEC3 p2   = (VEC3){ p2x, p2y, d_project_beam_to_sphere( TRACKBALLSIZE, p2x, p2y ) };
    VEC3 axis = METHOD(TPRE3, cross( p2, p1) );
    VEC3 d    = METHOD(TPRE3, sub(p1, p2)    );

    TYPE t = sqrt( 1 - METHOD(TPRE3, norm2(d) ) );
    
    return METHOD(TPREF, fromCosAngleAxis( t, axis ) );
}

//static inline VEC METHOD(TPREF, norm( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )
//static inline VEC METHOD(TPREF, cos ( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )
//static inline VEC METHOD(TPREF, dist( VEC v, TYPE x, TYPE y, TYPE z){ v.x*=x; v.y*=y; v.z*=z; return v; }   )

#undef  TPREF
#undef  TPRE3
#undef  TYPE
#undef  VEC
#undef  VEC3
#undef  MAT3