
#ifndef AeroSurf_h
#define AeroSurf_h

#include "fastmath.h"
#include "Vec3.h"

//#include <SDL2/SDL_opengl.h>
//#include <Draw3D.h>

#include "Body.h"

double const lowSpeedCuoff = 1e-8;

class AeroSurfaceDebugRecord { public:
    double ca;
    double sa;
    double CD;
    double CL;
    Vec3d  uair;
    Vec3d  force;
    Vec3d  gdpos;
};


struct AeroPolar {
	double CD0    = 0.02;
	double dCD    = 0.9;
	double dCDS   = 0.9;
	double dCL    = 6.28;
	double dCLS   = 2.82743338823;
	double sStall = 0.16;
	double wStall = 0.08;

    inline void getLD( double ca, double sa, double& CD, double& CL )const{
		//double abs_sa = fabs( sa );
		//printf("polarModel %f  %f %f %f  %f %f\n", area, CD0, dCD, dCDS, dCL, dCLS);
		double abs_sa = (sa>0)?sa:-sa;
		double wS     = trashold_cub<double>( abs_sa, sStall, sStall+wStall );
		double mS     = 1 - wS;
		CD     = CD0 + ( mS*dCD*abs_sa + wS*dCDS        ) * abs_sa;
		//CL     =       ( mS*dCL        + wS*dCLS*abs_ca ) * sa;
		if( ca <0 ){ ca=-ca; sa=-sa; };
		CL     =       ( mS*dCL        + wS*dCLS*ca ) * sa;
	}

};


class AeroSurface : public KinematicBody {
	public:
	RigidBody *craft;

	Vec3d C;     // Aerodynamic coefficients in each direction
	//Vec3d lpos;
	//Mat3d lrot;

    static constexpr double SAFETY_v = 1e-6;
	double area   = 1.0;  // [m^2]
	double CD0    = 0.02;
	double dCD    = 0.9;
	double dCDS   = 0.9;
	double dCL    = 6.28;
	double dCLS   = 2.82743338823;
	double sStall = 0.16;
	double wStall = 0.08;

	bool useC=false, usePolar=true;

    //bool DEBUGsurf = false;

    AeroSurfaceDebugRecord* dbgRec = NULL;


	// ==== function declarations

	void render( );
	void fromString( const char * str );
	void fromStringPolarModel( const char * str );

	// ==== inline functions

    inline void polarModel( double ca, double sa, double& CD, double& CL ){
		//double abs_sa = fabs( sa );
		//printf("polarModel %f  %f %f %f  %f %f\n", area, CD0, dCD, dCDS, dCL, dCLS);
		double abs_sa = (sa>0)?sa:-sa;
		double wS     = trashold_cub<double>( abs_sa, sStall, sStall+wStall );
		double mS     = 1 - wS;
		CD     = CD0 + ( mS*dCD*abs_sa + wS*dCDS        ) * abs_sa;
		//CL     =       ( mS*dCL        + wS*dCLS*abs_ca ) * sa;
		if( ca <0 ){ ca=-ca; sa=-sa; };
		CL     =       ( mS*dCL        + wS*dCLS*ca ) * sa;
	}

	inline Vec3d getPos(){ Vec3d gp; globalPos( craft->pos, craft->rotMat, gp ); return gp; };
	//inline Mat3d getRot(){ };

	inline void  applyForce( const Vec3d& vair0 ){
        //printf ( " %lf %lf %lf   %lf %lf %lf\n", lpos.x, lpos.y, lpos.z, C.a, C.b, C.c );

		Mat3d grot;
		//grot.set_mmul_NT( lrot, craft->rotMat );
		grot.set_mmul( lrot, craft->rotMat );
		//Mat3d rotMat; rotMat.setT(craft->rotMat);  grot.set_mmul( lrot, rotMat );
		Vec3d gdpos; craft->rotMat.dot_to_T( lpos, gdpos );

		//globalPos( {0.0,0.0} );
		//globalRot( {} );

		Vec3d uair;
		uair.set_cross( gdpos, craft->omega ); uair.add( vair0 );
		//uair.set( vair0 );
		//printf( "uair (%f,%f,%f) ", uair.x,uair.y,uair.z );

		double vrair2  = uair.norm2();
		if( vrair2 > lowSpeedCuoff ){
			double vrair = sqrt(vrair2);
			uair.mul( 1/vrair );

			// plane force
			double ca = grot.a.dot( uair );
			double cb = grot.b.dot( uair );
			double cc = grot.c.dot( uair );
			double prefactor = vrair2*area;

			Vec3d force;

			if( usePolar ){
                double CD,CL;
                polarModel( -cc, cb, CD, CL );
                //if(DEBUGsurf) printf("(%3.3f,%3.3f)  (%3.3f,%3.3f) \n", cc, cb, CD, CL );
                if(dbgRec){ dbgRec->ca=-cc; dbgRec->sa=cb; dbgRec->CD=CD; dbgRec->CL=CL; };
                CL*=prefactor; CD*=prefactor;

                double cb2 = cb*cb;
                if( (cb2>0.000001)&&(cb2<0.999999) ){
                    Vec3d airUp;
                    airUp.set_add_mul(grot.b,uair,-cb);
                    //printf(  " airUp (%g,%g,%g) cabc (%g,%g,%g) \n", airUp.x,airUp.y,airUp.z,  ca,cb,cc );
                    airUp.normalize();
                    force.set_lincomb( CL, airUp, CD, uair );

                    //printf(  " uair (%g,%g,%g) airUp (%g,%g,%g) \n", uair.x, uair.y, uair.z, airUp.x,airUp.y,airUp.z);
                    //printf(  " CLD %g,%g force (%g,%g,%g) \n", CL,CD, force.x,force.y,force.z);


                    //glColor3f( 1.0f,0.0f,0.0f ); Draw3D::drawVecInPos( uair   ,    craft->pos + gdpos );
                    //glColor3f( 0.0f,0.0f,1.0f ); Draw3D::drawVecInPos( grot.b ,    craft->pos + gdpos );
                    //glColor3f( 0.0f,1.0f,0.0f ); Draw3D::drawVecInPos( airUp  ,    craft->pos + gdpos );

                    //glColor3f( 0.0f,1.0f,0.0f ); Draw3D::drawVecInPos( airUp*CL  ,    craft->pos + gdpos );
                    //glColor3f( 0.0f,1.0f,0.0f ); Draw3D::drawVecInPos( airUp*CL  ,    craft->pos + gdpos );
                }else{
                    force.set_mul( uair, CD );
                }
			}else{
                force.set(0.0d);
			}

			if( useC ){
                force.add_mul( grot.a, C.a*ca*prefactor );
                force.add_mul( grot.b, C.b*cb*prefactor );
                force.add_mul( grot.c, C.c*cc*prefactor );
            }
			//printf( "vrair %f \n", vrair );
			//printVec( uair ); printf("uair\n");

			if(dbgRec){ dbgRec->force=force;  dbgRec->gdpos=gdpos; dbgRec->uair=uair*sqrt(vrair2); };
			//printf( "force (%f,%f,%f) \n", force.x,force.y,force.z );
			craft->apply_force( force, gdpos );

			//drawMatInPos( craft->rotMat, craft->pos );
			//Draw3D::drawMatInPos( grot, craft->pos + gdpos );

			//glColor3f( 0.9f,0.0f,0.9f ); Draw3D::drawVecInPos( force     *0.005,    craft->pos + gdpos );
			//glColor3f( 0.0f,0.5f,0.5f ); Draw3D::drawVecInPos( uair*vrair*0.1,   craft->pos + gdpos );

		}
	};

	//inline void setAoA( const Mat3d& lrot0, double AoA ){};

	// TODO: this can be probably deleted ?
	/*
	inline void  applyForceSimple( const Vec3d& vair0 ){
        //printf ( " %lf %lf %lf   %lf %lf %lf\n", lpos.x, lpos.y, lpos.z, C.a, C.b, C.c );

		Mat3d grot;
		//Mat3d rotMat; rotMat.setT(craft->rotMat);  grot.set_mmul( lrot, rotMat );
		grot.set_mmul_NT( lrot, craft->rotMat );
		//grot.set_mmul( lrot, craft->rotMat );
		//grot.set_mmul( craft->rotMat, lrot );

		Vec3d gdpos; craft->rotMat.dot_to( lpos, gdpos );

		//printf( "gdpos    %3.3f %3.3f %3.3f \n", gdpos.x, gdpos.y, gdpos.z );
		//printf( "grot.a %3.3f %3.3f %3.3f \n", grot.ax, grot.ay, grot.az );
		//printf( "grot.b %3.3f %3.3f %3.3f \n", grot.bx, grot.by, grot.bz );
		//printf( "grot.c %3.3f %3.3f %3.3f \n", grot.cx, grot.cy, grot.cz );

        //printf( "lrot.a %3.3f %3.3f %3.3f \n", lrot.ax, lrot.ay, lrot.az );
		//printf( "lrot.b %3.3f %3.3f %3.3f \n", lrot.bx, lrot.by, lrot.bz );
		//printf( "lrot.c %3.3f %3.3f %3.3f \n", lrot.cx, lrot.cy, lrot.cz );



        //glColor3f( 1.0f,0.0f,0.0f ); Draw3D::drawVecInPos( grot.a*C.a,    craft->pos + gdpos );
        //glColor3f( 0.0f,1.0f,0.0f ); Draw3D::drawVecInPos( grot.b*C.b,    craft->pos + gdpos );
        //glColor3f( 0.0f,0.0f,1.0f ); Draw3D::drawVecInPos( grot.c*C.c,    craft->pos + gdpos );

		Vec3d uair;
		uair.set_cross( gdpos, craft->omega );
		//uair.set(0);
		uair.add( vair0 );
		double vrair2  = uair.norm2();
		if( vrair2 >0 ){
			double vrair = sqrt(vrair2);
			uair.mul( 1/vrair );

			// plane force
			double ca = grot.a.dot( uair );
			double cb = grot.b.dot( uair );
			double cc = grot.c.dot( uair );
			double prefactor = vrair2*area;

			Vec3d force;
			//force.set( uair*(C.y*vrair2) );
			//force.set( grot.b*(C.y*cb*vrair2) );
			force.set( grot.a*(C.a*ca*prefactor) + grot.b*(C.b*cb*prefactor) + grot.c*(C.c*cc*prefactor) );

			//printf( "vrair %f \n", vrair );
			//printVec( uair ); printf("uair\n");

			craft->apply_force( force, gdpos );

			//drawMatInPos( craft->rotMat, craft->pos );
			//Draw3D::drawMatInPos( grot, craft->pos + gdpos );

			//glColor3f( 0.9f,0.0f,0.9f ); Draw3D::drawVecInPos( force     *0.005,    craft->pos + gdpos );
			//glColor3f( 0.0f,0.5f,0.5f ); Draw3D::drawVecInPos( uair*vrair*0.1,   craft->pos + gdpos );

		}
	};
	*/

};

#endif






