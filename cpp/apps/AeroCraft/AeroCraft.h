
#ifndef AeroCraft_h
#define AeroCraft_h

#include "fastmath.h"
#include "Body.h"
#include "AeroSurf.h"


#include <SDL2/SDL_opengl.h>
#include <Draw3D.h>


class Propeler{
    public:
    Vec3d lpos;
    Vec3d dir;
    double area;
    double power;
    double efficiency;
    double CD;
    double vstatic; // velocity of stream for properler stationary with respect to air

    inline double getThrust( double v0 ){
        // derived from those equations
        // dm = S * v0
        // F = dm * Dv
        // P = F * v0 + 0.5*dm*(Dv**2) = S*(v0**2)*Dv + 0.5*S*v0*(Dv**2)
        double dm = area*(v0+vstatic);
        double a  = 0.5*dm;
        double b  = dm*v0;
        double c  = -power;
        double Dv1,Dv2;
        quadratic_roots( a, b, c,  Dv1, Dv2 );
        return dm * Dv1 * efficiency - dm*v0*CD;
    }

    inline void fromString( const char * str ){
        sscanf ( str, " %lf %lf %lf    %lf %lf %lf    %lf %lf %lf %lf",
                    &lpos.x, &lpos.y, &lpos.z,
                    &dir.x, &dir.y, &dir.z,
                    &area, &power, &efficiency, &CD
                );
        dir.normalize();
        vstatic = pow(4*power/area,0.333333);
        printf ( "%lf %lf %lf %lf\n", area, power, efficiency, CD );
    }

};


class AeroCraft : public RigidBody {
	public:

	int nPanels = 0, nPropelers= 0;
	AeroSurface * panels    = NULL;
	Propeler    * propelers = NULL;

	AeroSurface * leftAirelon =NULL;
	AeroSurface * rightAirelon=NULL;
	AeroSurface * elevator    =NULL;
	AeroSurface * rudder      =NULL;

	Vec3d totalThrust;

	double maxAileron  = 0.1;
	double maxElevator = 0.5;
	double maxRudder   = 0.5;

	// ==== function declarations

	virtual void render();

	int fromFile( const char * fname );

	// ==== inline functions

	inline void applyAeroForces( const Vec3d& vwind ){
		Vec3d vair = vwind - vel;

        //panels[0].DEBUGsurf = true;
		for( int i=0; i<nPanels; i++ ){
            //printf( " %i %i \n", i, nPanels );
            //panels[i].applyForceSimple( vair );
            panels[i].applyForce( vair );
		}

		//Mat3 rmat; rmat.setT(rotMat);
		totalThrust.set(0.0d);
		for( int i=0; i<nPropelers; i++ ){
            Vec3d gdpos,gdir;
            rotMat.dot_to( propelers[i].dir,  gdir  );
            rotMat.dot_to( propelers[i].lpos, gdpos );
            double vdir   = -gdir.dot( vair );  // printf("vdir %g \n", vdir);
            double thrust = propelers[i].getThrust(vdir);
            //glColor3d(1.0f,1.0f,0.0f); Draw3D::drawVecInPos( gdir, gdpos+pos );  Draw3D::drawPointCross( gdpos+pos, 0.5 );
            gdir.mul(thrust);
            //glColor3d(1.0f,1.0f,0.0f); Draw3D::drawVecInPos( gdir, gdpos+pos );  Draw3D::drawPointCross( gdpos+pos, 0.2 );
            totalThrust.add(gdir);
            apply_force( gdir, gdpos );
		}
	}


	inline void steerTo(double droll, double dpitch, double dyaw){
        panels[0].lrot.rotate(  _clamp( -droll ,-maxAileron ,maxAileron ), {1.0,0.0,0.0} );
        panels[1].lrot.rotate(  _clamp( +droll ,-maxAileron ,maxAileron ), {1.0,0.0,0.0} );
        panels[2].lrot.rotate(  _clamp(  dpitch,-maxElevator,maxElevator), {1.0,0.0,0.0} );
        panels[3].lrot.rotate(  _clamp( -dyaw  ,-maxRudder  ,maxRudder  ), {0.0,1.0,0.0} );
    }

};



class AeroCraftControler{
    public:
    AeroCraft * craft;

    double vvert_target    =-10.0d;
    //double vvert_rate      =0.05d;
    //double vvert_strength  =0.000005d;
    //double vvert_strength  =0.00005d;  // this is safe without "vy_against_fy" checking
    double vvert_strength  =0.005d;
    double dvvert_smooth   =0.0d;

	inline void control( double dt ){

        double dvvert  = (vvert_target - craft->vel.y);
        //double bmix    = vvert_rate*dt;
        //dvvert_smooth  = bmix*dvvert + (1-bmix)*dvvert_smooth;
        dvvert_smooth = dvvert*dt;

        bool vy_against_fy = (( craft->force.y < 0 )&&( dvvert_smooth < 0 ))  && !(( craft->force.y > 0 )&&( dvvert_smooth > 0 ));
        if( !vy_against_fy ){
            //printf("%g %g %g %g\n", craft->vel.y, dvvert, bmix, dvvert_smooth );
            craft->elevator->lrot.rotate(  _clamp(  dvvert_smooth*vvert_strength,-craft->maxElevator,craft->maxElevator), {1.0,0.0,0.0} );
        }

	}

};

#endif
