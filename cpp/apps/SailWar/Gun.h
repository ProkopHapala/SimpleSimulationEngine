
#ifndef Gun_h
#define Gun_h

#include "Vec3.h"
#include "Body.h"

class Projectile;

class Gun : public KinematicBody {
	public:
	double muzzle_velocity = 1.0;
	double projectile_mass = 1.0;

	void set_direction( Vec3d dir );
	Projectile * fireProjectile( const Vec3d& gpos, const Mat3d& grot, const Vec3d& gvel  );

};

#endif 

