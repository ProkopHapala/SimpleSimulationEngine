
#ifndef AeroSurf_h
#define AeroSurf_h

#include "fastmath.h"
#include "Vec3.h"

#include <SDL2/SDL_opengl.h>
#include <drawMath.h> 

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
		Mat3d grot;  grot.set_mmul( lrot, craft->rotMat );
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
			force.set( grot.a*(C.x*ca*vrair2) + grot.b*(C.y*cb*vrair2) + grot.c*(C.z*cc*vrair2) );

			//printf( "vrair %f \n", vrair );
			//printVec( uair ); printf("uair\n");
			//printVec( uair ); printf("force\n");

			craft->apply_force( force, gdpos );

			//drawMatInPos( craft->rotMat, craft->pos );
			drawMatInPos( grot, craft->pos + gdpos );

			//glColor3f( 0.0f,0.0f,0.0f );
			//drawLine( craft->pos + gdpos,craft->pos + gdpos + (grot.b*5));
			glColor3f( 0.9f,0.0f,0.9f );
			drawLine( craft->pos + gdpos,craft->pos + gdpos + (force*0.1));
			//glColor3f( 0.0f,0.5f,0.0f );
			//drawLine( craft->pos + gdpos,craft->pos + gdpos + (uair*2));

		}
	};

};

#endif






