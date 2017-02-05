
#ifndef Battleship_h
#define Battleship_h

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>

#include "Vec3.h"
#include "Ship2D.h"
#include "Projectile3D.h"
#include "Gun3D.h"
#include "Turret.h"
#include "Collisions.h"

class Battleship : public Ship2D, public CollisionObject {
	public:
    char * name;
	int nguns;
	Turret * turrets;

	double life_max               = 1.0d;
	double life                   = 1.0d;
	double life_regeneration_rate = 2.9d;

	double reload_rate   = 2.5d;
	double gunload_left  = 1.0d;
	double gunload_right = 1.0d;

	// ==== function declarations

	Turret ** initGuns( int n, Vec3d pos1, Vec3d pos2, Vec3d ldir, double muzzle_velocity );
	//void initAllGuns( int n );
	//void drawGun( Gun * gun );
	//virtual void draw( );
	//virtual void drawHitBox( );
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, Vec3d * where, Vec3d * normal );
	virtual bool loadFromFile( char const* filename );
    virtual void update( double dt );

};

#endif  // #ifndef Battleship_h

