#ifndef Projectile3D_h
#define Projectile3D_h

#include <vector>

#include "fastmath.h"
#include "Vec3.h"
#include "Body.h"
#include "geom3D.h"

#include "ShooterCommon.h"

class Projectile3D : public PointBody { public:
	int    id,kind;
    Vec3d  old_pos; // DEPRECATED: use pos+vel*dt does the job
    double time=0;

    inline void update_Projectile3D( double dt ){
        time   += dt;
        old_pos.set(pos);
        move_PointBody(dt);
    }

    virtual void update( double dt ){
        update_Projectile3D( dt );
    };

    // TO DO : would be usefull if material of target specified
    //virtual void hit(){}

};

class Burst3d { public:
    int    id=-1;
    double age=0;
    bool   discard=false;    // NOT NECESSARY  ... we may simply shift it to     1e+300 or to NaN

    ProjectileType* type = 0;

    Capsula3D bbox;
    Box       aabb; // We may have easily two bounding boxes

	std::vector<Particle3d> shots;


    void updateBBox( double dt ){
        Vec3d hdir;
        bbox.p = shots[0].pos;
        //hdir = shots[0].vel; hdir.normalize(); // STRATEGY 1
        Vec3d op0; shots.back().getOldPos(dt,op0); hdir.set_sub(shots[0].pos,op0); hdir.normalize(); // STRATEGY 2
        double lmin = +1e+300;
        double lmax = -1e+300;
        float r2max=0.0;
        int n = shots.size();
        //printf( "p0 (%g,%g,%g) p1 (%g,%g,%g) \n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z );
        for( int i=0; i<n; i++ ){
            Vec3d d;
            double r2,l;
            //Vec3d p = shots[i].pos;
            const Particle3d& p = shots[i];
            d.set_sub( p.pos-p.vel*-dt, bbox.p); l=d.makeOrthoU(hdir); r2=d.norm2(); if(r2>r2max) r2max=r2; if(l>lmax)lmax=l; if(l<lmin)lmin=l; // printf(" fw %i %g %g \n", i, l, lmin);  //d.add_mul( hdir, -hdir.dot(d) );
            p.getOldPos(dt,d);    d.sub(bbox.p); l=d.makeOrthoU(hdir); r2=d.norm2(); if(r2>r2max) r2max=r2; if(l>lmax)lmax=l; if(l<lmin)lmin=l;
        }
        bbox.r = sqrt(r2max);
        bbox.l = lmax-lmin;
        //printf( "lmin %g lmax %g \n", lmin, lmax );
        bbox.p.add_mul( hdir, lmin );
        bbox.hdir=hdir;
	}


	void updateBBox_torq( double dt ){
        // Explanation
        //  - We can evaluate how much tilted are particles with respect hdir-axis by averaging 1st-moment
        //    Sum_i{(d_i-d_av)*(l_i-l_av)} = Sum_i{ d_i*l_i - d_av*l_i - d_i*l_av + d_av*l_av }  = Sum_i {  d_i*l_i - d_av*l_av }
        //  - The problem is that after corresponding rotation we need to recalculate r,l bounds => we would do it only in next step
        Vec3d hdir;
        bbox.p = shots[0].pos;
        //hdir = shots[0].vel; hdir.normalize(); // STRATEGY 1
        //Vec3d op0; shots.back().getOldPos(dt,op0); hdir.set_sub(shots[0].pos,op0); hdir.normalize(); // STRATEGY 2
        double lmin = +1e+300;
        double lmax = -1e+300;
        float r2max=0.0;
        int n = shots.size();
        //printf( "p0 (%g,%g,%g) p1 (%g,%g,%g) \n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z );
        Vec3d torq = Vec3dZero, dav = Vec3dZero;
        double lav =0;
        for( int i=0; i<n; i++ ){
            Vec3d d;
            double r2,l;
            //Vec3d p = shots[i].pos;
            const Particle3d& p = shots[i];
            d.set_sub( p.pos-p.vel*-dt, bbox.p); l=d.makeOrthoU(hdir); torq.add_mul(d,l); dav.add(d); lav+=l; r2=d.norm2(); if(r2>r2max) r2max=r2; if(l>lmax)lmax=l; if(l<lmin)lmin=l; // printf(" fw %i %g %g \n", i, l, lmin);  //d.add_mul( hdir, -hdir.dot(d) );
            p.getOldPos(dt,d);    d.sub(bbox.p); l=d.makeOrthoU(hdir); torq.add_mul(d,l); dav.add(d); lav+=l; r2=d.norm2(); if(r2>r2max) r2max=r2; if(l>lmax)lmax=l; if(l<lmin)lmin=l; // printf(" bk %i %g %g \n", i, l, lmin);
        }
        double invN = 1.0/shots.size();
        torq.mul(invN);
        torq.add_mul( dav, -lav*invN*invN );
        torq.mul(invN);

        bbox.r = sqrt(r2max);
        bbox.l = lmax-lmin;
        //printf( "lmin %g lmax %g \n", lmin, lmax );
        bbox.p.add_mul( hdir, lmin );
        bbox.hdir=hdir;
	}

    void move(double dt, const Vec3d& accel0, double airDensity ){
        //printf("Burst3d::move \n");
        if(discard) return;
        int i=0;
        int n = shots.size();
        //Vec3d tmpPos[n];
        double balisticCoef = airDensity * type->balisticCoef;
        //printf( "balisticCoef %g \n", balisticCoef );
        // TODO : we can calculate p.vel.norm() and supersonic drag by taylor expansion
        // double rv0 = shots[0].vel.norm();
        for( int i=0; i<n; i++ ){

            Particle3d& p = shots[i];
            //printf( "shot[%i] p (%g,%g,%g) v (%g,%g,%g)\n", i, p.pos.x, p.pos.y, p.pos.z, p.vel.x, p.vel.y, p.vel.z );

            //if( (p.vel.norm()<100.0) ) printf( "shot[%i] p (%g,%g,%g) v (%g,%g,%g)\n", i, p.pos.x, p.pos.y, p.pos.z, p.vel.x, p.vel.y, p.vel.z );

            //tmpPos[i] = p.pos;
            Vec3d accel;
            //accel=accel0;
            accel.set_add_mul( accel0, p.vel, p.vel.norm()*-balisticCoef );
            p.move(dt, accel );
            //p.move(dt, accel0 );
            age+=dt;
        }
        updateBBox( dt ); // we asume shot[0] is most forward, shot[n-1] is least

        //printf( "BBox %g %g %g \n", bbox.p.x, bbox.p.y, bbox.p.z );
        //printf("Burst3d::move DONE\n");
    }

    void addShot( const Vec3d& pos, const Vec3d& vel ){
        shots.push_back( (Particle3d){0.0,pos,vel} );
        //Particle3d& p = shots.back(); printf( "add shot[%i] p (%g,%g,%g) v (%g,%g,%g)\n", shots.size(), p.pos.x, p.pos.y, p.pos.z, p.vel.x, p.vel.y, p.vel.z );
    }

    inline void hit( int i){
        shots[i].age = 1e+300;
    };

    Burst3d(){};
    Burst3d( ProjectileType* type_, int id_ ){ type=type_; id=id_; }

};

#endif
