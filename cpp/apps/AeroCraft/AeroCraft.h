
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

	// ==== function declarations

	virtual void render();

	// ==== inline functions

	inline void applyAeroForces( const Vec3d& vwind ){
		Vec3d vair = vwind - vel;
		wingLeft .applyForceSimple( vair );
		wingRight.applyForceSimple( vair );
		rudder   .applyForceSimple( vair );
		elevator .applyForceSimple( vair );
	};

};

#endif
