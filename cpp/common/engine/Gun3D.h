
#ifndef Gun3D_h
#define Gun3D_h

#include "Vec3.h"
#include "Body.h"

#include "Projectile3D.h"

class Projectile3D;

/*
class Gun3d{ public:
    Vec3d dir;
    double muzzle_velocity = 1.0;
    ProjectileType* prjType = 0;

    Burst3d* burst = 0;

};
*/

class Gun3D : public KinematicBody { public:
	double muzzle_velocity = 1.0;
	double Projectile3D_mass = 1.0;

	void set_direction( Vec3d dir );
	Projectile3D * fireProjectile3D( const Vec3d& gpos, const Mat3d& grot, const Vec3d& gvel  );

};

#endif

