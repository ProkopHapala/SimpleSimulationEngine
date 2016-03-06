
#ifndef Yacht2D_h
#define Yacht2D_h

#include "Vec2.h"
#include "Body2D.h"
#include "AeroSurf2D.h"

class SailWarWorld;

class Yacht2D : public RigidBody2D {
	public:
	SailWarWorld * world;
	AeroSurf2D keel;
	AeroSurf2D rudder;
	AeroSurf2D mast;

	// ==== function declarations

	virtual void draw( );
	virtual void applySailForces( const Vec2d& windSpeed, const Vec2d& watterSpeed );

    virtual void testSail    ( int n, int nsub, double dt, double windSpeed, Vec2d * poss, Vec2d * vels, Vec2d * rots );
    virtual int  convergeSail( int nmax, int nsub, double dt, double windSpeed, double vtol, double rtol, Vec2d& vel_conv, Vec2d& rot_conv );
    virtual void evalPolar   ( int n, double dt, double vtol, double rtol, double * phi_rudder, double * phi_mast, double * wind_speed,   Vec2d * vels, Vec2d * rots, bool reseting );

};

#endif  // #ifndef Yacht2D_h

