
#ifndef SailWarWorld_h
#define SailWarWorld_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "Convex2d.h"

#include "Projectile.h"
#include "Projectile.h"
#include "Frigate2D.h"

class SailWarWorld {
	public:
	double ground_level;
	Vec3d  wind_speed;
	Vec2d  watter_speed;

	int perFrame = 10;
	double dt = 0.0001;

	int defaultShipShape;
	CollisionShape * defaultCollisionShape;

	std::vector<Convex2d*>   isles;

	std::vector<Frigate2D*>  ships;
	std::vector<Projectile*> projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector

    void makeShip( const Vec2d& pos, double angle, char * filename, int shape, CollisionShape * collisionShape );
	void update_world ( );
	void init_world   ( );

//	void projectile_collisions();

};

#endif  // #ifndef SailWarWorld_h
