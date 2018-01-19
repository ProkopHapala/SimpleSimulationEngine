
#ifndef GridIndex2D_h
#define GridIndex2D_h

#include "fastmath.h"
#include "Vec2.h"

class GridIndex2D{ public:
    Vec2i     n    = (Vec2i){0,0};
    int       ntot = 0;
    inline void   setN(Vec2i n_)       { n=n_; ntot=n.x*n.y;                }
    inline Vec2i  i2ip(int   i ) const { return {i%n.x,i/n.x};   } // https://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder
    inline int    ip2i(Vec2i ip) const { return (n.x*ip.y+ip.x); }
    inline bool   validIndex( Vec2i ip)const{ return (ip.x>0)&&(ip.x<n.x)&&(ip.y>0)&&(ip.y<n.y); }
    inline Vec2i  wrap_index( Vec2i ip)const{ return (Vec2i){ wrap_index_fast(ip.a,n.a),  wrap_index_fast(ip.b,n.b) }; }
    inline Vec2i  clam_index( Vec2i ip)const{ return (Vec2i){clamp_index_fast(ip.a,n.a), clamp_index_fast(ip.b,n.b) }; }
    inline double fetchWraped( Vec2i ip, const double * hs)const{ int ia=wrap_index_fast(ip.a,n.a); int ib=wrap_index_fast(ip.b,n.b); return hs[n.x*ip.y+ip.x]; }
};

#endif

