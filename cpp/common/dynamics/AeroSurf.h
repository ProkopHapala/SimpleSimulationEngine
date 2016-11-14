
#ifndef AeroSurf_h
#define AeroSurf_h

#include "fastmath.h"
#include "Vec3.h"

#include <SDL2/SDL_opengl.h>
#include <Draw3D.h>

#include "Body.h"

class AeroSurface : public KinematicBody {
	public:
	RigidBody *craft;

	Vec3d C;     // Aerodynamic coefficients in each direction
	Vec3d lpos;
	Mat3d lrot;


    static constexpr double SAFETY_v = 1e-6;
	double area;  // [m^2]
	double CD0    = 0.02;
	double dCD    = 0.9;
	double dCDS   = 0.9;
	double dCL    = 6.28;
	double dCLS   = 2.82743338823;
	double sStall = 0.16;
	double wStall = 0.08;

	// ==== function declarations

	void render( );
	void fromString( const char * str );

	// ==== inline functions

    inline void polarModel( double ca, double sa, double& CD, double& CL ){
		//double abs_sa = fabs( sa );
		double abs_sa = (sa>0)?sa:-sa;
		double wS     = trashold_cub<double>( abs_sa, sStall, sStall+wStall );
		double mS     = 1 - wS;
		CD     = CD0 + ( mS*dCD*abs_sa + wS*dCDS        ) * abs_sa;
		//CL     =       ( mS*dCL        + wS*dCLS*abs_ca ) * sa;
		if( ca <0 ){ ca=-ca; sa=-sa; };
		CL     =       ( mS*dCL        + wS*dCLS*ca ) * sa;
	}

	inline void  applyForceSimple( const Vec3d& vair0 ){
        //printf ( " %lf %lf %lf   %lf %lf %lf\n", lpos.x, lpos.y, lpos.z, C.a, C.b, C.c );

		Mat3d grot; Mat3d rotMat;
		rotMat.setT(craft->rotMat);  grot.set_mmul( lrot, rotMat );
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

			Vec3d force;
			//force.set( uair*(C.y*vrair2) );
			//force.set( grot.b*(C.y*cb*vrair2) );
			force.set( grot.a*(C.a*ca*vrair2) + grot.b*(C.b*cb*vrair2) + grot.c*(C.c*cc*vrair2) );

			//printf( "vrair %f \n", vrair );
			//printVec( uair ); printf("uair\n");

			craft->apply_force( force, gdpos );

			//drawMatInPos( craft->rotMat, craft->pos );
			//Draw3D::drawMatInPos( grot, craft->pos + gdpos );




			glColor3f( 0.9f,0.0f,0.9f ); Draw3D::drawVecInPos( force     *0.005,    craft->pos + gdpos );
			glColor3f( 0.0f,0.5f,0.5f ); Draw3D::drawVecInPos( uair*vrair*0.1,   craft->pos + gdpos );

		}
	};






	inline void  applyAeroForce( const Vec3d& vair0 ){
        //printf ( " %lf %lf %lf   %lf %lf %lf\n", lpos.x, lpos.y, lpos.z, C.a, C.b, C.c );

		Mat3d grot; Mat3d rotMat;
		rotMat.setT(craft->rotMat);  grot.set_mmul( lrot, rotMat );
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
		if( vrair2 < 1e-16 ) return;

        double vrair = sqrt(vrair2);
        uair.mul( 1/vrair );

        double sa = uair.dot( rotMat.b );
        double ca = uair.dot( rotMat.c );
        double CD,CL;
        polarModel( ca, sa, CD, CL );
        // TO DO :   polar model works only if  uair is close to parallel with rotMat.c
        // this means that ca~1.0 ... we can use this as an interpolation parameter,
        // and for the rest stereo-angle use simple drag of plane

        Vec3d ulift,uslide, force;
        uslide.set_cross(uair,rotMat.b); uslide.normalize();
        ulift .set_cross(uair,uslide);   ulift .normalize();


        force.set_lincomb( CD*area*vrair2, CL*area*vrair2,  0,     uair, ulift, uslide );

        craft->apply_force( force, gdpos );

        //drawMatInPos( craft->rotMat, craft->pos );
        //Draw3D::drawMatInPos( grot, craft->pos + gdpos );

        glColor3f( 0.9f,0.0f,0.9f ); Draw3D::drawVecInPos( force     *0.005,    craft->pos + gdpos );
        glColor3f( 0.0f,0.5f,0.5f ); Draw3D::drawVecInPos( uair*vrair*0.1,   craft->pos + gdpos );

	};


};

#endif






