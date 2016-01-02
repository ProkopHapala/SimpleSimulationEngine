
#ifndef Frigate2D_h
#define Frigate2D_h

//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Vec3.h"
#include "Yacht2D.h"
#include "Projectile.h"
#include "Gun.h"
#include "Collisions.h"

class Frigate2D : public Yacht2D, public CollisionObject {
	public:
    char * name;
	int nguns;
	Gun ** left_guns;
	Gun ** right_guns;

	double life_max               = 1.0d;
	double life                   = 1.0d;
	double life_regeneration_rate = 2.9d;

	double reload_rate   = 2.5d;
	double gunload_left  = 1.0d;
	double gunload_right = 1.0d;

	// ==== function declarations

	Gun ** initGuns( int n, Vec3d pos1, Vec3d pos2, Vec3d ldir, double muzzle_velocity );
	void initAllGuns( int n );
	void fire_gun_row( int n, Gun ** guns, std::vector<Projectile*> * projectiles );
	void fire_left ( std::vector<Projectile*> * projectiles );
	void fire_right( std::vector<Projectile*> * projectiles );
	void drawGun( Gun * gun );
	virtual void draw( );
	virtual void drawHitBox( );
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, Vec3d * where, Vec3d * normal );
	bool loadFromFile( char const* filename );
    virtual void update( double dt );

};

#endif  // #ifndef Frigate2D_h

