
#ifndef SoftPolyLine3D_h
#define SoftPolyLine3D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
//#include "geom3D.h"


class SoftPolyLine3d{ public:
    //std::vector<PointBody2D> ps;
    int n;
    //Quat4d * qrots = NULL;
    Vec3d  * ps    = NULL;
    Vec3d  * vs    = NULL;
    Vec3d  * fs    = NULL;

    // aux
    Vec3d  * ups   = NULL;
    Vec3d  * fws   = NULL;
    double * ls    = NULL;

    double*  mass = NULL; // mass
    double*  l0s  = NULL;
    double * kls  = NULL; // stiffness in lenght
    double * kas  = NULL; // stiffness in angle
    double * kts  = NULL; // stiffness in torq

    void allocate( int n_, bool tempF, bool tempH, bool bParams ){
        n  = n_;
        ps = new Vec3d[n];
        vs = new Vec3d[n];
        if(tempF) fs = new Vec3d[n];
        if(tempH){
            //printf("aloc tempH\n");
            fws = new Vec3d [n-1];
            ups = new Vec3d [n-1];
            ls  = new double[n-1];
        }
        if(bParams){
            mass  = new double[n];
            l0s   = new double[n-1];
            kls   = new double[n-1];
            kas   = new double[n-2];
            kts   = new double[n-2];
        }
    }

/*
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
*/

    void cleanForce   (          ){ for(int i=0; i<n; i++){ fs[i]=Vec3dZero; } }
    void mulVelocities( double f ){ for(int i=0; i<n; i++){ vs[i].mul(f);    } }

    // TODO: move this is to radial forces (?)
    void updateAux(){
        Vec3d op=ps[0];
        for(int i=0; i<n-1; i++){
            Vec3d p  = ps[i+1];
            Vec3d h  = p-op;
            ls [i]   = h.normalize();
            fws[i]   = h;
            Vec3d up = ups[i];
            up.add_mul( h, -up.dot(h) );
            up.normalize();
            ups[i] = up;
            op=p;
        }
    }

    void radialForces(){
        Vec3d of = Vec3dZero;
        for(int i=0; i<n-1; i++){
            Vec3d f; f.set_mul( fws[i] , ( l0s[i] - ls[i] ) * kls[i] );
            //Vec2d f; f.set_mul( hs[i], ( -(l0s[i] - ls[i]) ) );
            of   .sub(f);
            fs[i].add(of);
            //printf( "%i  (%f,%f) (%f,%f) %f \n", i, hs[i].x, hs[i].y, fs[i].x, fs[i].y, f );
            of=f;
        }
        fs[n-1].add(of);
    }

    void angularForces(){
        Vec3d oup = ups[0];
        Vec3d of  = Vec3dZero;
        Vec3d oof = Vec3dZero;
        double ol = ls[0];
        for(int i=0; i<n-2; i++){
            Vec3d& up = ups[i+1];
            Vec3d& fw = fws[i+1];
            double  l = ls [i+1];
            // calc torq
            Vec3d  tq;   tq .set_cross( up, oup );
            double cfw = fw.dot(tq);
            tq.add_mul( fw, -cfw      );
            tq.mul    (         kts[i] );
            tq.add_mul( fw, cfw*kas[i] );
            // calc force
            Vec3d fa,fb;
            fa.set_cross( oup, tq*(1/ol) );
            fb.set_cross(  up, tq*(1/ l) );
            // store force
            fs[i] = oof+fa;
            oof.sub(fa); oof.sub(fb);
            of .add(fb);
        }
        fs[n-2].add(oof);
        fs[n-1].add(of);
    }

/*
// ---- chain with fixed lengths
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
*/

/*
// ---- distance from point
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
*/

    void getForceNoTemp(double dt){}

    void move(double dt, double damp){
        for(int i=0; i<n; i++){
            Vec3d& v = vs[i];
            v.mul(1-damp);
            v.add_mul( fs[i], dt/mass[i] );
            ps[i].add_mul( v, dt );
        }
    }

    int nearestPoint( Vec3d p ){
        int imin=0; double r2min=1e+300;
        for(int i=0; i<n; i++){
            double r2=(ps[i]-p).norm2();
            if(r2<r2min){r2min=r2; imin=i;}
        }
        return imin;
    }

    void placePointsLine( Vec3d p0, Vec3d p1 ){
        Vec3d dp = (p1-p0)*(1.0/(n-1));
        for(int i=0; i<n; i++){
            ps  [i] = p0 + dp*i;
            printf( "%i (%f,%f) \n", i, ps[i] );
        }
    }

    // set current gometry as equlibrium geom
    void setCurrentRef(){
        Vec3d op=ps[0];
        Vec3d oh;
        for(int i=1; i<n; i++){
            Vec3d p  = ps[i];
            Vec3d h  = p-op;
            l0s  [i-1] = h.normalize();
            //if(i>1)ang0s[i-2].set_udiv_cmplx(oh,h);
            // TODO - angular reference

            oh=h; op=p;
            //printf("%i a0: (%f,%f) l: %f\n", i, ang0s[i-2].x, ang0s[i-2].y, l0s[i-1] );
        }
    }

    void init(int n_, Vec3d p0, Vec3d p1, double m, double kl, double ka ){
        allocate( n_, true, true, true );
        //printf("DEBUG 1 \n");
        placePointsLine( p0, p1 );
        //printf("DEBUG 2 \n");
        for(int i=0; i<n  ; i++){ mass[i] = m; vs[i].set(0.0); } //printf("DEBUG 3 \n");
        for(int i=0; i<n-1; i++){ kls [i] = kl; }// printf("DEBUG 3 \n");
        for(int i=0; i<n-2; i++){ kas [i] = ka; } //printf("DEBUG 4 \n");
        setCurrentRef(); //printf("DEBUG 5 \n");
    }

};

#endif

