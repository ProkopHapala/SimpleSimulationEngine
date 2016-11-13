
#ifndef AeroCraft_h
#define AeroCraft_h

#include "Body.h"
#include "AeroSurf.h"

class AeroCraft : public RigidBody {
	public:
	AeroSurface wingLeft;
	AeroSurface wingRight;
	AeroSurface elevator;
	AeroSurface rudder;

	int nPanels = 0;
	AeroSurface * panels = NULL;

	// ==== function declarations

	virtual void render();

	int fromFile( const char * fname );

	// ==== inline functions

	inline void applyAeroForces( const Vec3d& vwind ){
		Vec3d vair = vwind - vel;
		/*
		wingLeft .applyForceSimple( vair );
		wingRight.applyForceSimple( vair );
		rudder   .applyForceSimple( vair );
		elevator .applyForceSimple( vair );
        */
		for( int i=0; i<nPanels; i++ ){
            //printf( " %i %i \n", i, nPanels );
            panels[i].applyForceSimple( vair );
		}
	};

};

#endif
