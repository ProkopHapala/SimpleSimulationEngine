
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
    OMesh *  mesh = NULL;
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
        vphi        = vphi*DEG2RAD;
        vtheta      = vtheta*DEG2RAD;

        OMesh * mesh_ = new OMesh();
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

	double phi_target, theta_target;

	double reload      = 1.0;
    Vec3d  gun_dir;

	//Vec3d target;

    virtual void fromString( char const* str ){
        sscanf( str, "%i %lf %lf %lf %lf %lf %lf %lf\n", &kind,  &lpos.x, &lpos.y, &lpos.z,    &phimin, &phimax,  &thetamin, &thetamax );
        printf(      "%i %lf %lf %lf %lf %lf %lf %lf\n",  kind,   lpos.x,  lpos.y,  lpos.z,     phimin,  phimax,   thetamin,  thetamax );
        thetamin *= DEG2RAD;  thetamax *= DEG2RAD;
        phimin   *= DEG2RAD;  phimax   *= DEG2RAD;
        phi   = 0.5*(phimin+phimax);
        theta = 0.5*(thetamin+thetamax);
    }

	void update(double dt ){
        if(reload<1) reload += dt*type->reload_rate;
        phi   += clamp_abs( dangle( phi_target-phi ),   dt*type->vphi   );
        theta += clamp_abs( theta_target-theta, dt*type->vtheta );

        //printf( "theta_target %f theta %f \n", theta_target, theta );
	};

    virtual void updateTransforms( const Vec3d& pos0, const Mat3d& rot0 ){
        //Mat3d lrot;
        //lrot.fromEuler( 0.0, 0.0, 0.0 );
        //lrot.print(); exit(0);
        //lrot.setOne();
        //lrot.rotate_csa( sin(phi), cos(phi), {0.0,1.0,0.0} );
        //grot.set_mmul( lrot, rot0 );
        rot0.dot_to( lpos, gpos ); gpos.add( pos0 );
        //phi=1.0;
        grot.set(rot0);

        //grot.print();
        //grot.setOne();
        grot.rotate_csa( cos(phi), sin(phi), {0.0,1.0,0.0} );

        grot.dot_to( {cos(theta), sin(theta), 0.0 }, gun_dir );
	}


	void aim( const Vec3d& target, double G, double phi0 ){
        Vec3d d; d.set_sub(target,gpos);
        double dist    = sqrt(d.x*d.x+d.z*d.z);
        //phi_target     = atan2(d.x,d.z) + M_PI_2 + phi0;
        phi_target     = atan2(d.x,d.z) + M_PI_2 + M_PI + phi0;
        //phi_target     = atan2(d.x,d.z) + M_PI_2;

        // https://en.wikipedia.org/wiki/Trajectory_of_a_projectile    #Angle θ {\displaystyle \theta } \theta required to hit coordinate (x,y)
        double x1,x2;
        double a = -0.5*(-G)*sq(dist/type->vmuzzle);
        quadratic_roots( a, dist, -d.y+a, x1, x2 );
        if(x2<x1)x1=x2;
        theta_target = atan( x1 );

        //printf( "  x1 %f x2 %f theta_target %f \n", x1, x2, theta_target );
	}



	//void   set_direction( Vec3d dir );
	//Projectile3D * fireProjectile3D( const Vec3d& gvel );
	int shoot( std::vector<Projectile3D*>& projectiles, const Vec3d& gvel );

	//virtual void draw();

};

#endif

