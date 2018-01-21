
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


class SoftPolyLine2D{ public:
    //std::vector<PointBody2D> ps;
    int n;
    Vec2d * ps = NULL;
    Vec2d * vs = NULL;
    Vec2d * fs = NULL;

    Vec2d  * hs = NULL; // unitary bond direction vectors
    double * ls = NULL; // bond direction vectors

    double * mass  = NULL; // mass
    double * l0s   = NULL; // relaxed lengths
    double * kls   = NULL; // stiffnes length
    Vec2d  * ang0s = NULL; // zero rotation
    double * kas   = NULL; // stiffness angular

    void allocate( int n_, bool tempF, bool tempH, bool bParams ){
        n  = n_;
        ps = new Vec2d[n];
        vs = new Vec2d[n];
        if(tempF) fs = new Vec2d[n];
        if(tempH){
            //printf("aloc tempH\n");
            hs = new Vec2d [n-1];
            ls = new double[n-1];
        }
        if(bParams){
            //printf("aloc bParams\n");
            mass  = new double[n]; // mass
            l0s   = new double[n-1]; // relaxed lengths
            kls   = new double[n-1]; // stiffnes length
            ang0s = new Vec2d [n-2]; // zero rotation
            kas   = new double[n-2]; // stiffness angular
        }
    }

    //void initFF( double l0, double  ){}

    //#define lh(p,op, l,h)   h=p-op; l=h.normalize();
    //#define

    void setCurrentRef(){
        Vec2d op=ps[0];
        Vec2d oh;
        for(int i=1; i<n; i++){
            Vec2d p  = ps[i];
            Vec2d h  = p-op;
            l0s  [i-1] = h.normalize();
            if(i>1)ang0s[i-2].set_udiv_cmplx(oh,h);
            oh=h; op=p;
            //printf("%i a0: (%f,%f) l: %f\n", i, ang0s[i-2].x, ang0s[i-2].y, l0s[i-1] );
        }
    }

    void cleanForce   (          ){ for(int i=0; i<n; i++){ fs[i].set(0.0); } }
    void mulVelocities( double f ){ for(int i=0; i<n; i++){ vs[i].mul(f);   } }

    // TODO: move this is to radial forces (?)
    void updateAux(){
        Vec2d op=ps[0];
        for(int i=0; i<n-1; i++){
            Vec2d p = ps[i+1];
            Vec2d h = p-op;
            ls[i]   = h.normalize();
            hs[i]   = h;
            //printf( "%i   (%f,%f) %f \n", i, hs[i].x, hs[i].y, ls[i] );
            //Vec2d::mul_cmplx(op,d);
            op=p;
        }
    }

    void radialForces(){
        Vec2d of; of.set(0.0);
        for(int i=0; i<n-1; i++){
            Vec2d f; f.set_mul( hs[i] , ( l0s[i] - ls[i] ) * kls[i] );
            //Vec2d f; f.set_mul( hs[i], ( -(l0s[i] - ls[i]) ) );
            of   .sub(f);
            fs[i].add(of);
            //printf( "%i  (%f,%f) (%f,%f) %f \n", i, hs[i].x, hs[i].y, fs[i].x, fs[i].y, f );
            of=f;
        }
        fs[n-1].add(of);
    }

    void angularForces(){
        Vec2d oh; oh.set_perp(hs[0]);
        Vec2d of;  of.set(0.0);
        Vec2d oof; oof.set(0.0);
        double ol = ls[0];
        for(int i=1; i<n-1; i++){
            Vec2d h,h_,fa,fb;
            h.set_perp( hs[i] );
            h_.set_mul_cmplx(h,ang0s[i-1]);
            double torq = h_.cross(oh) * kas[i-1];
            double l=ls[i];
            fa.set_mul(oh,torq/ol); // f[i-1]
            fb.set_mul(h ,torq/l ); // f[i+1]
            //fs[i-1].add(fa);
            //fs[i  ].sub(fa+fb);
            //fs[i+1].add(fb);
            //oof.add(fa);
            of     .sub(fa+fb);
            fs[i-1].add(oof+fa);
            oof=of; of=fb;
            //glColor3f(0.0,0.0,1.0); Draw2D::drawVecInPos_d( oh, ps[i] );
            //glColor3f(1.0,0.0,0.0); Draw2D::drawVecInPos_d( h , ps[i] );
            //glColor3f(0.0,0.5,0.0); Draw2D::drawVecInPos_d( h_, ps[i] );
            oh=h; ol=l;
            //of=f;
        }
        fs[n-2].add(oof);
        fs[n-1].add(of);
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
        for(int i=1; i<n; i++){
            Vec2d  h  = hs[i-1];
            double c  = h.dot(odp);
            Vec2d  dp = p-ps[i];
            //if(c>0){glColor3f(0.0,0.0,1.0);}else{glColor3f(1.0,0.0,0.0);}
            //Draw2D::drawLine_d( p, ps[i-1] );
            //glColor3f(0.0,1.0,0.0); Draw2D::drawVecInPos_d( h  *100.0, ps[i-1] );
            //glColor3f(1.0,0.0,1.0); Draw2D::drawVecInPos_d( odp*10.0, ps[i-1] );
            if( c>0 ){
                if( c<ls[i-1] ){ // dist from line segement
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
        return dpmin;
    }

    void getForceNoTemp(double dt){}

    void move(double dt, double damp){
        for(int i=0; i<n; i++){
            Vec2d v = vs[i];
            Vec2d f = fs[i];
            //double cvf    = v.dot(f);
            //double fdamp = 1.0;
            //if(damp>0){ fdamp -= damp*cvf/sqrt( v.norm2() * f.norm2() ); } // velocity dependent force
            //v.add_mul( f, fdamp * dt/mass[i] );
            v.mul(1-damp);
            v.add_mul( f, dt/mass[i] );
            vs[i]=v;
            ps[i].add_mul( v, dt );
        }
    }

    int nearestPoint( Vec2d p ){
        int imin=0; double r2min=1e+300;
        for(int i=0; i<n; i++){
            double r2=(ps[i]-p).norm2();
            if(r2<r2min){r2min=r2; imin=i;}
        }
        return imin;
    }

    /*
    void update(double dt){
        if( hs ){
            updateAux();

        }else{ moveNoTemp(dt) }
    }
    */

    void placePoints( Vec2d p0, Vec2d p1 ){
        Vec2d dp = (p1-p0)*(1.0/(n-1));
        for(int i=0; i<n; i++){
            ps  [i] = p0 + dp*i;
            printf( "%i (%f,%f) \n", i, ps[i] );
        }
    }

    void init(int n_, Vec2d p0, Vec2d p1, double m, double kl, double ka ){
        allocate( n_, true, true, true );
        //printf("DEBUG 1 \n");
        placePoints( p0, p1 );
        //printf("DEBUG 2 \n");
        for(int i=0; i<n  ; i++){ mass[i] = m; vs[i].set(0.0); } //printf("DEBUG 3 \n");
        for(int i=0; i<n-1; i++){ kls [i] = kl; }// printf("DEBUG 3 \n");
        for(int i=0; i<n-2; i++){ kas [i] = ka; } //printf("DEBUG 4 \n");
        setCurrentRef(); //printf("DEBUG 5 \n");
    }

};

#endif

