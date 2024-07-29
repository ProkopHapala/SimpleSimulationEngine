
#ifndef  Vec3Utils_h
#define  Vec3Utils_h

#include "Vec3.h"
#include "Mat3.h"


template<class T> static inline void evalLenghs( int n, Vec2i* bonds, const Vec3T<T>* ps, T* lengths ){ for(int i=0; i<n; i++){ Vec2i b=bonds[i]; lengths[i]=ps[b.i].dist(ps[b.j]); } }

template<class T> 
static inline T getCOG( int n, const Vec3T<T>* vs, const T* ws, Vec3T<T>& cog ){
    cog.set(0.0);
    T wsum=0;
    for(int i=0; i<n; i++){ T w=ws[i]; wsum+=w; cog.add_mul(vs[i],w); }
    cog.mul( 1/wsum );
    return wsum;
}

//namespace Vec3{
template<class T> 
static inline Vec3T<T> sum( int n, const Vec3T<T>* vs ){
    Vec3T<T> v; v.set(0.0);
    for(int i=0; i<n; i++){ v.add(vs[i]); }
    return v;
}

template<class T> 
static inline Vec3T<T> sum( int n, const Vec3T<T>* vs, int* selection ){
    Vec3T<T> v; v.set(0.0);
    for(int i=0; i<n; i++){ v.add(vs[selection[i]]); }
    return v;
}
template<class T> static inline Vec3T<T> average( int n, const Vec3T<T>* vs                 ){ Vec3T<T> v=sum(n,vs          ); v.mul( 1/(T)n ); return v; }
template<class T> static inline Vec3T<T> average( int n, const Vec3T<T>* vs, int* selection ){ Vec3T<T> v=sum(n,vs,selection); v.mul( 1/(T)n ); return v; }

template<class T>
static inline double findNearest( int n, const Vec3T<T>* ps, const Vec3d& p0, int& ifound, int ino=-1 ){
    double r2min=1e+300; int imin=-1;
    for(int i=0; i<n; i++){ if(i==ino)continue; Vec3d d; d.set_sub(ps[i],p0); double r2=d.norm2(); if(r2<r2min){r2min=r2;imin=i;} }
    ifound=imin;
    return r2min;
}

template<class T> 
static inline T maxR2( int n, const Vec3T<T>* vs ){
    T R2=0;
    for(int i=0; i<n; i++){ double r2=vs[i].norm2(); if(r2>R2)R2=r2; }
    return R2;
}


template<class T>
static inline Vec3T<T> torq( int n, const Vec3T<T>* ps,const Vec3T<T>* vs, const Vec3T<T>& p0, int* selection=0 ){
    Vec3T<T> v; v.set(0.0);
    for(int i=0; i<n; i++){ 
        int ii=i; if(selection)ii=selection[i];
        Vec3T<T> dp; dp.set_sub( ps[ii], p0 );
        v.add_cross( dp, vs[ii] );
    }
    return v;
}
template<class T>
static inline Vec3T<T> torq( int n, const Vec3T<T>* ps, const Vec3T<T>* vs ){
    Vec3T<T> p0 = average(n,ps );
    return torq(n, ps, vs, p0 );
}

template<class T>
static inline void apply_torq( int n, Vec3T<T> p0, Vec3T<T> ax, const Vec3T<T>* ps, const Vec3T<T>* vs ){
    for(int i=0; i<n; i++){
        Vec3d dp = ps[i] - p0;
        Vec3d v; v.set_cross(ax,dp);
        vs[i].add(v);
    }
}

template<class T> inline void move  (int n, Vec3T<T>* ps, const Vec3T<T>& shift                                     ){ for(int i=0; i<n; i++)ps[i].add(shift); }
template<class T> inline void scale (int n, Vec3T<T>* ps, const Vec3T<T>& center, const Vec3T<T>& sc                ){ for(int i=0; i<n; i++){ ps[i].scale     ( sc, center      ); };  }
template<class T> inline void rotate(int n, Vec3T<T>* ps, const Vec3T<T>& center, const Vec3T<T>& uaxis, T ca, T sa ){ for(int i=0; i<n; i++){ ps[i].rotate_csa( ca, sa, uaxis, center   ); }; }
template<class T> inline void rotate(int n, Vec3T<T>* ps, const Vec3T<T>& center, const Vec3T<T>&  axis, T phi      ){ rotate( n,ps,center,axis.normalized(), cos(phi), sin(phi) ); }

template<class T> inline void move  (int n, const int* selection, Vec3T<T>* ps, const Vec3T<T>& shift                                     ){ for(int i =0;  i<n;  i++)  ps[selection[i]].add(shift);                        }
template<class T> inline void scale (int n, const int* selection, Vec3T<T>* ps, const Vec3T<T>& center, const Vec3T<T>& sc                ){ for(int ii=0; ii<n; ii++){ ps[selection[ii]].scale     ( sc, center ); }; }
template<class T> inline void rotate(int n, const int* selection, Vec3T<T>* ps,       Vec3T<T>  center, const Vec3T<T>& uaxis, T ca, T sa ){ for(int ii=0; ii<n; ii++){ ps[selection[ii]].rotate_csa( ca, sa, uaxis, center   ); }; }
template<class T> inline void rotate(int n, const int* selection, Vec3T<T>* ps, const Vec3T<T>& center, const Vec3T<T>&  axis, T phi      ){ rotate( n,selection,ps,center,axis.normalized(), cos(phi), sin(phi) ); }

template<class T>
static inline void affineTransform( int n, const Vec3T<T>* vin, Vec3T<T>* vout, const Mat3T<T>& oM, const Mat3T<T>& newM, const Vec3T<T> oP=Vec3T<T>{0,0,0}, const Vec3T<T> newP=Vec3T<T>{0,0,0} ){
    Mat3d M1,M; 
    oM.invert_to( M1 );
    M.set_mmul(newM,M1);
    for(int i=0; i<n; i++){ 
        Vec3d v; M.dot_to( vin[i]-oP, v ); v.add(newP); vout[i]=v;
    }
}

template<class T>
void orient( int n, Vec3T<T>* ps, const Vec3T<T>& p0, const Vec3T<T>& dir, const Vec3T<T>& up ){
    Mat3T<T> rot; rot.fromDirUp( dir, up );
    for(int i=0; i<n; i++){
        Vec3T<T> p = ps[i]-p0;
        rot.dot_to(p,ps[i]);
        ps[i].add(p0);
    }
}

template<class T> 
static inline void bbox( Vec3T<T>& pmin, Vec3T<T>& pmax, int n, const Vec3T<T>* ps, const int* selection=0, bool bInit=true ){
    if(bInit){
        pmin=(Vec3T<T>){+1e+37,+1e+37,+1e+37};
        pmax=(Vec3T<T>){-1e+37,-1e+37,-1e+37};
    }
    for(int i=0; i<n; i++){ 
        int ii=i;
        if(selection)ii=selection[i];
        pmin.setIfLower  (ps[ii]);
        pmax.setIfGreater(ps[ii]);
    }
}

template<class T> 
static inline Vec3T<T> cog_bbox( int n, const Vec3T<T>* ps, const int* selection=0 ){
    Vec3d pmin,pmax;
    bbox( pmin, pmax, n, ps, selection, true );
    return (pmin+pmax)*0.5;
}

template<class T> 
void cogTo0( int n, Vec3T<T>* ps, const int* selection=0 ){
    Vec3T<T>& pmin,pmax;
    bbox(pmin,pmax, n, ps, selection, true );
    pmin.add(pmax); pmin.mul(0.5);
    move(n,ps,pmin);
}

template<class T> 
void findRotation( Vec3T<T>& p0, int n, const Vec3T<T>* ps, Mat3T<T>& rot){
    Mat3T<T> XX=(Mat3T<T>){0,0,0};
    for(int i=0; i<n; i++){
        Vec3T<T> dp; dp.set_sub(ps[i],p0);
        XX.addOuter( dp, dp, 1.0 );
    }
    //printf( "XX: " ); printMat(XX);
    Vec3T<T> evs;
    XX.eigenvals(evs);
    evs.sort();
    //printf(  "FindRotation evs(%g,%g,%g) \n", evs.x,evs.y,evs.z );
    Vec3T<T> c;
    XX.eigenvec( evs.x, rot.a );
    XX.eigenvec( evs.y, rot.b );
    XX.eigenvec( evs.z, rot.c );
}

template<typename T, typename Func>
void numDeriv( Vec3T<T> p, T d, Vec3T<T>& f, Func func){
    T d_=d*0.5;
    p.x+=d_; f.x = func(p); p.x-=d; f.x-=func(p); p.x+=d_;
    p.y+=d_; f.y = func(p); p.y-=d; f.y-=func(p); p.y+=d_;
    p.z+=d_; f.z = func(p); p.z-=d; f.z-=func(p); p.z+=d_;
    f.mul(1/d);
}

template<typename T>
void makeSamples(const Vec2i& ns, const Vec3T<T>& p0, const Vec3T<T>& a, const Vec3T<T>& b, Vec3T<T> *ps ){
    Vec3T<T> da=a*(1.0/ns.x);
    Vec3T<T> db=b*(1.0/ns.y);
    //printf( "da (%g,%g,%g)\n", da.x,da.y,da.z );
    //printf( "db (%g,%g,%g)\n", db.x,db.y,db.z );
    for(int ib=0; ib<ns.y; ib++){
        Vec3T<T> p = p0+db*ib;
        for(int ia=0; ia<ns.x; ia++){
            *ps = p;
            p.add(da);
            ps++;
        }
    }
}

//}

#endif