
#ifndef Turret_h
#define Turret_h

#include "Projectile3D.h"
#include "Gun3D.h"
#include "Object3D.h"

#include "Vec3.h"
#include "Body.h"

//class Projectile3D;

class Turret : public Object3D {
	public:
	double muzzle_velocity   = 1.0;
	double Projectile3D_mass = 1.0;
	double phimin  =-1000.0,phimax  =1000.0;
	double thetamin=-5.0,   thetamax=80.0;

	//void   set_direction( Vec3d dir );
	Projectile3D * fireProjectile3D( const Vec3d& gpos, const Mat3d& grot, const Vec3d& gvel  );

};

#endif

