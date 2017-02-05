
#ifndef Ship2D_h
#define Ship2D_h

#include "Vec2.h"
#include "Body2D.h"
#include "AeroSurf2D.h"
#include "Mesh.h"
#include "Warrior25D.h"

//class GameWorld;

class Ship2D : public Warrior25D {
	public:
	//GameWorld * world;
	AeroSurf2D keel;
	AeroSurf2D rudder;
	AeroSurf2D mast;
	//int Sails;
	//AeroSurf2D * Sail;      // TO DO
	//Propeler propeler;      // TO DO

	// ==== function declarations

	//virtual void loadFromFile( char const* fname ) = 0;
	//virtual void fromString( char const* str );

	virtual void update( double dt, const Vec3d& wind_speed  );

	void applyHydroForces( const Vec2d& watterSpeed );
	void applySailForces ( const Vec2d& windSpeedd  );
	void updatePropeler  ( const Vec2d& watterSpeed );

    virtual void testSail    ( int n, int nsub, double dt, double windSpeed, Vec2d * poss, Vec2d * vels, Vec2d * rots );
    virtual int  convergeSail( int nmax, int nsub, double dt, double windSpeed, double vtol, double rtol, Vec2d& vel_conv, Vec2d& rot_conv );
    virtual void evalPolar   ( int n, double dt, double vtol, double rtol, double * phi_rudder, double * phi_mast, double * wind_speed,   Vec2d * vels, Vec2d * rots, bool reseting );

    inline double getPropelerThrust( double v0 ){
        /*
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
        */
        return power * throttle / sqrt( vel.norm2() + v0*v0 ) ;
    }
};

#endif  // #ifndef Ship2D_h

