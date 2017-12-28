
#ifndef AeroCraft_h
#define AeroCraft_h

#include "fastmath.h"
#include "Body.h"
#include "AeroSurf.h"

//#include <SDL2/SDL_opengl.h>
//#include <Draw3D.h>

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
	AeroSurface * panels      =NULL;
	Propeler    * propelers   =NULL;

	AeroSurface * leftAirelon =NULL;
	AeroSurface * rightAirelon=NULL;
	AeroSurface * elevator    =NULL;
	AeroSurface * rudder      =NULL;

	Vec3d totalThrust;

	double maxAileron  = 0.1;
	double maxElevator = 0.5;
	double maxRudder   = 0.5;

	// ==== function declarations

	//virtual void render();

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

class SurfControl{ public:

    double val       = 0.0;
    double maxRate   = 0.01;
    double relaxRate = maxRate*0.2;
    double minval    =-0.1;
    double maxval    = 0.1;

    AeroSurface * surf;

    void setSymetricRange( double mval ){ minval=-mval; maxval=mval; }
    void apply  ( double dval ){ surf->lrot.rotate( dval, surf->lrot.a );  val+=dval;  }
    void inc  (){ double dval=maxval-val; dval=_min(dval, maxRate); apply( dval ); }
    void dec  (){ double dval=minval-val; dval=_max(dval,-maxRate); apply( dval ); }
    void relax(){ double dval=_clamp(-val,-relaxRate,relaxRate);       apply( dval ); }

};

class AeroCraftControler{
    public:
    AeroCraft * craft  = NULL;
    AeroCraft * craft0 = NULL;

    double vvert_target    = 0.0d;
    double vvert_strength  = 0.1d;
    double dvvert_smooth   = 0.0d;

    double roll_target     = 0.9d;
    //double roll_strength   = 0.1d;
    //double roll_damp       = 50;
    double roll_strength   = 0.3d;
    double roll_damp       = 30;

    SurfControl leftAirelon;
    SurfControl rightAirelon;
    SurfControl rudder;
    SurfControl elevator;

    void attach( AeroCraft * craft_ ){
        craft             =craft_;
        leftAirelon .surf =craft->leftAirelon;
        rightAirelon.surf =craft->rightAirelon;
        rudder      .surf =craft->rudder;
        elevator    .surf =craft->elevator;
    }


	inline void control( double dt ){
        double dvvert  = (vvert_target - craft->vel.y);
        bool vy_against_fy = (craft->force.y * dvvert) < 0.0d;
        if( vy_against_fy ){
            //double bmix    = vvert_rate*dt;
            //dvvert_smooth  = bmix*dvvert + (1-bmix)*dvvert_smooth;
            dvvert_smooth = dvvert;
            double dpitch = dvvert_smooth*vvert_strength*dt;
            //printf("%f %f\n", dpitch, dt );
            //printf("%g %g %g %g\n", craft->vel.y, dvvert, bmix, dvvert_smooth );
            craft->elevator->lrot.rotate(  _clamp( dpitch,-craft->maxElevator,craft->maxElevator), {1.0,0.0,0.0} );
        }
        //
        Mat3d rmat; rmat.setT(craft->rotMat);
        double roll_err = (roll_target - rmat.a.y);
        double comega   =  craft->omega.dot( rmat.c );
        bool adjust_roll = true;
        //if( fabs(comega>1e-8) )  adjust_roll = ((roll_err * comega ) < 0.0);
        roll_err *= (1- roll_damp*roll_err*comega);
        double droll=0;
        if( adjust_roll ){
            droll = _clamp( roll_err*roll_strength*dt,-craft->maxAileron,craft->maxAileron);
            craft-> leftAirelon->lrot.rotate(  droll , {1.0,0.0,0.0} );
            craft->rightAirelon->lrot.rotate( -droll , {1.0,0.0,0.0} );
        }
        //printf("autoPilot roll_err %g droll %g    %g \n",roll_err, droll,  roll_err * comega );
	}

    void resetSteer( ){
        for( int i=0; i<craft->nPanels; i++ ){ craft->panels[i].lrot=craft0->panels[i].lrot; }
    }

    void steerToDir( const Vec3d& dir ){
        Mat3d rotMatT;
        rotMatT.setT(craft->rotMat);
        //Draw3D::drawMatInPos( rotMatT, craft->pos );
        double dyaw   = rotMatT.a.dot( dir );
        double dpitch = rotMatT.b.dot( dir );
        const double acut = 0.1;
        //double droll = (a>acut)?(a-acut):((a<-acut)?(a+acut):0.0d);
        double droll = dyaw;
        resetSteer();
        craft->steerTo( 0.1*droll, 0.5*dpitch+0.5*fabs(dyaw), 0.5*dyaw);
    };

};

#endif
