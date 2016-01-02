
#ifndef GameWorld_h
#define GameWorld_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"

#include "Projectile.h"
#include "Projectile.h"
#include "Frigate2D.h"

class GameWorld {
	public:
	double ground_level;
	Vec3d  wind_speed;
	Vec2d  watter_speed;

	int perFrame = 10;
	double dt = 0.0001;

	std::vector<Frigate2D*>  ships;
	std::vector<Projectile*> projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector

	void update( );
	void init( );

//	void projectile_collisions();

};

#endif  // #ifndef GameWorld_h
