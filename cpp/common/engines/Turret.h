
#ifndef Turret_h
#define Turret_h

#include "Projectile3D.h"
#include "Gun3D.h"
#include "Object3D.h"

#include "Vec2.h"
#include "Vec3.h"
#include "Body.h"

//class Projectile3D;


class TurretType{
    public:
    Mesh *  mesh = NULL;
    int     shape   = 0;
    int     ngun    = 1;
    double  prjmass = 1.0;
    double  vmuzzle = 1000.0;
    double  reload_rate = 0.1;
    double  vphi    = 10.0;   // deg/s
    double  vtheta  = 10.0;   // deg/s
    double  mass    = 10e+3;  //
    double  Rsize   = 10.0;


    void fromString( char const* str ){
        char meshname[256];
        sscanf( str, "%s %lf %lf %i %lf %lf %lf %lf %lf\n", meshname, &Rsize, &mass, &ngun,  &prjmass, &vmuzzle, &reload_rate, &vphi, &vtheta );

        printf(      "%s %lf %lf %i %lf %lf %lf %lf %lf\n", meshname,  Rsize,  mass,  ngun,   prjmass,  vmuzzle,  reload_rate,  vphi,  vtheta );

        reload_rate = 1/reload_rate;
        vphi        = M_PI*vphi/180.0;
        vtheta      = M_PI*vtheta/180.0;

        Mesh * mesh_ = new Mesh();
        int res = mesh_->fromFileOBJ( meshname );
        if(res>0) mesh=mesh_;
    }

};


class Turret : public Object3D {
	public:
	//double muzzle_velocity   = 1.0;
	//double Projectile3D_mass = 1.0;
	TurretType * type = NULL;
	double phi,theta;
	double phimin  =-1000.0,phimax  =1000.0;
	double thetamin=-5.0,   thetamax=80.0;

	double reload      = 1.0;
    Vec3d  gun_dir;

	//Vec3d target;

    virtual void fromString( char const* str ){
        sscanf( str, "%i %lf %lf %lf %lf %lf %lf %lf\n", &kind,  &lpos.x, &lpos.y, &lpos.z,    &phimin, &phimax,  &thetamin, &thetamax );
        printf(      "%i %lf %lf %lf %lf %lf %lf %lf\n",  kind,   lpos.x,  lpos.y,  lpos.z,     phimin,  phimax,   thetamin,  thetamax );
        phi   = 0.5*(phimin+phimax);
        theta = 0.5*(thetamin+thetamax);
    }

	void update(double dt ){
        if(reload<1) reload += dt*type->reload_rate;
	};

    virtual void updateTransforms( const Vec3d& pos0, const Mat3d& rot0 ){
        //Mat3d lrot;
        lrot.fromEuler( 0, phi, 0 );

        rot0.dot_to( lpos, gpos ); gpos.add( pos0 );
        grot.set_mmul( lrot, rot0 );

        grot.dot_to( {cos(theta), sin(theta)}, gun_dir );
	}

	//void   set_direction( Vec3d dir );
	Projectile3D * fireProjectile3D( const Vec3d& gvel );

	//virtual void draw();

};

#endif

