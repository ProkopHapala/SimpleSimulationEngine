
#ifndef Projectile_h
#define Projectile_h

#include "Vec3.h"
#include "Body.h"

class GameWorld;

class Projectile : public PointBody {
	public:
	GameWorld * world;
	double dragCoef;
	double penetration;
	double damage;
	Vec3d  old_pos;

	void ground_hit( );
	bool check_hit(  );
	void addDragForce( const Vec3d& vwind, Vec3d& aeroForce );
	virtual void evalForce( );
	virtual void move( double dt );
	virtual void draw();
};

#endif 
