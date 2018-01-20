
#ifndef SoftPolyLine2D_h
#define SoftPolyLine2D_h

#include "fastmath.h"
#include "Vec2.h"

/*

TODO:

INSPIRATION:
   - simulation of projectile impact on 2d plate. Several different plates
*/

/*
template <typename T> struct Vec2{T x,y;};
template <typename T> T cross(Vec2<T> a,Vec2<T> b){ return a.x*b.y - a.y*b.x; }
template <typename T> T dot  (Vec2<T> a,Vec2<T> b){ return a.x*b.y + a.y*b.x; }
*/

/*
template <typename T>
struct Vec2{
    using V = Vec2<T>;
    T x,y;
    // functions
    static T cross(V a,V b){ return a.x*b.y - a.y*b.x; }
    static T dot  (V a,V b){ return a.x*b.y + a.y*b.x; }
};

using Vec2dd = Vec2<double>;
*/




/*
template <typename T>{
    struct Vec2{T x,y;};
    T cross(Vec2<T> a,Vec2<T> b){ return a.x*b.y - a.y*b.x; }
    T dot  (Vec2<T> a,Vec2<T> b){ return a.x*b.y + a.y*b.x; }
}
*/



/*
template <typename T> T mul  (Vec2<T> a,Vec2<T> b){ return (Vec2<T>){a.x*b.y,a.y*b.x; }
template <typename T> T add  (Vec2<T> a,Vec2<T> b){ return (Vec2<T>){a.x+b.y,a.y+b.x; }
template <typename T> T sub  (Vec2<T> a,Vec2<T> b){ return (Vec2<T>){a.x+b.y,a.y+b.x; }
template <typename T> T div  (Vec2<T> a,Vec2<T> b){ return (Vec2<T>){a.x/b.y,a.y/b.x; }
*/


class SoftPolyLine2D{
    //std::vector<PointBody2D> ps;
    int n;
    Vec2d * ps = NULL;
    Vec2d * vs = NULL;
    Vec2d * fs = NULL;

    Vec2d  * hs = NULL; // unitary bond direction vectors
    double * ls = NULL; // bond direction vectors

    double * l0s   = NULL; // relaxed lengths
    double * mass  = NULL; // mass
    double * ks    = NULL; // stiffness
    Vec2d  * ang0s = NULL; // zero rotation

    void allocate( int n, bool tempF, bool tempH ){
        n  = n;
        ps = new Vec2d[n];
        vs = new Vec2d[n];
        if(tempF) fs = new Vec2d[n];
        if(tempH){
            hs = new Vec2d [n-1];
            ls = new double[n-1];
        }
    }

    //void initFF( double l0, double  ){}

    void updateAux(){
        Vec2d op=ps[0];
        for(int i=1; i<n; i++){
            Vec2d p = ps[i];
            Vec2d d = p-op;
            ls[i-1] = d.normalize();
            hs[i-1] = d;
            //Vec2d::mul_cmplx(op,d);
        }
    }

    void radialForces(){
        Vec2d of;
        for(int i=0; i<n-1; i++){
            Vec2d f; f.set_mul( hs[i] , ( l0s[i] - ls[i] ) * ks[i] );
            of.sub(f);
            fs[i].add(of);
            of=f;
        }
        fs[n-1]=of;
    }

    void angularForces(){
        Vec2d oh; oh.set_perp(hs[0]);
        Vec2d of; of.set(0.0);
        for(int i=1; i<n-1; i++){
            Vec2d h,h_,f;
            h.set_perp( hs[i] );
            h_.set_mul_cmplx(h,ang0s[0]);
            double ff = h_.cross(oh) * ks[i];
            f .set_mul(h ,ff);
            of.add_mul(oh,ff);
            fs[i]=of;
            oh=h; of=f;
        }
        fs[n-1]=of;
    }

    void radialConstrain( int ipivot ){
        if( ipivot<(n-1) ){
            Vec2d op=ps[ipivot];
            for(int i=ipivot+1; i<n; i++){
                Vec2d  p=ps[i];
                Vec2d dp=p-op;
                double l = dp.norm();
                p.add_mul( dp, (l0s[i]-l)/l );
                ps[i]=p;
                op   =p;
            }
        }
        if( ipivot>0 ){
            Vec2d op=ps[ipivot];
            for(int i=ipivot+1; i>=0; i--){
                Vec2d  p=ps[i];
                Vec2d dp=p-op;
                double l = dp.norm();
                p.add_mul( dp, (l0s[i]-l)/l );
                ps[i]=p;
                op   =p;
            }
        }
    }

    Vec2d dpmin( const Vec2d& p ){ // derivative of distance squared
        //  line from point distance:
        //   ( (x2-x1)*y - (y2-y1)*x + x2*y1 - x1*y2 ) / |p1-p2|
        Vec2d  odp   = p-ps[0];
        Vec2d  dpmin = odp;
        double r2min = odp.norm2(); // dist first corner point
        for(int i=0; i<n-1; i++){
            Vec2d  dp = hs[i+1];
            Vec2d  h  = hs[i  ];
            double c  = h.dot(odp);
            if( c>0 ){
                if( c<l0s[i  ] ){ // dist from line segement
                    odp.add_mul(h,-c);
                    double r2 = odp.norm2();
                    if( r2<r2min ){ dpmin=odp; r2min=r2; }
                }else{             // dist from next corner point
                    double r2 = dp.norm2();
                    if( r2<r2min ){ dpmin=dp; r2min=r2; }
                }
            }
            odp=dp;
        }
    }

    void getForceNo(double dt){}

    void move(double dt){
        for(int i=0; i<n; i++){
            Vec2d v = vs[i];
            v.add_mul( fs[i], dt/mass[i] );
            vs[i]=v;
            ps[i].add_mul( v, dt );
        }
    }

    /*
    void update(double dt){
        if( hs ){
            updateAux();

        }else{ moveNoTemp(dt) }
    }
    */

};

#endif

