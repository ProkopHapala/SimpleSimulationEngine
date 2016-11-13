
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

	// ==== function declarations

	void render( );

	// ==== inline functions

	inline void  applyForceSimple( const Vec3d& vair0 ){
		Mat3d grot; Mat3d rotMat;
		rotMat.setT(craft->rotMat);  grot.set_mmul( lrot, rotMat );
		//grot.set_mmul( lrot, craft->rotMat );
		//grot.set_mmul( craft->rotMat, lrot );

		Vec3d gdpos; craft->rotMat.dot_to( lpos, gdpos );

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

			glColor3f( 1.0f,0.0f,0.0f ); Draw3D::drawVecInPos( grot.a*C.a,    craft->pos + gdpos );
			glColor3f( 0.0f,1.0f,0.0f ); Draw3D::drawVecInPos( grot.b*C.b,    craft->pos + gdpos );
			glColor3f( 0.0f,0.0f,1.0f ); Draw3D::drawVecInPos( grot.c*C.c,    craft->pos + gdpos );


			glColor3f( 0.9f,0.0f,0.9f ); Draw3D::drawVecInPos( force*0.1,    craft->pos + gdpos );
			glColor3f( 0.0f,0.5f,0.5f ); Draw3D::drawVecInPos( uair*vrair,   craft->pos + gdpos );

		}
	};

};

#endif






