
#ifndef Yacht2D_h
#define Yacht2D_h

#include "Vec2.h"
#include "Body2D.h"
#include "AeroSurf2D.h"

class Yacht2D : public RigidBody2D {
	public:
	AeroSurf2D keel;
	AeroSurf2D rudder;
	AeroSurf2D mast;

	// ==== function declarations
	
	virtual void draw( );
	virtual void applySailForces( const Vec2d& windSpeed, const Vec2d& watterSpeed );

};

#endif  // #ifndef Yacht2D_h

