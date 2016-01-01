
#ifndef Frigate2D_h
#define Frigate2D_h

//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Vec3.h"
#include "Yacht2D.h"
#include "Projectile.h"
#include "Gun.h"

class Frigate2D : public Yacht2D {
	public:
	
	int nguns;
	Gun ** left_guns; 
	Gun ** right_guns; 

	// ==== function declarations

	Gun ** initGuns( int n, Vec3d pos1, Vec3d pos2, Vec3d ldir, double muzzle_velocity );
	void initAllGuns( int n );
	void fire_gun_row( int n, Gun ** guns, std::vector<Projectile*> * projectiles );
	void fire_left ( std::vector<Projectile*> * projectiles );
	void fire_right( std::vector<Projectile*> * projectiles );
	void drawGun( Gun * gun );
	virtual void draw( );
	bool loadFromFile( char const* filename );

};

#endif  // #ifndef Frigate2D_h

