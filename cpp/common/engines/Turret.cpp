
#include "Projectile3D.h"
#include "Gun3D.h"

#include "Turret.h" // THE HEADER

/*
void Turret::set_direction( Vec3d dir ){
	Vec3d up; up.set( 0.0d, 0.0d, 1.0d );
	bounds.orientation.a.set( dir );
	bounds.orientation.c.set_cross( lrot.a, bounds.orientation.b  );
	bounds.orientation.b.set_cross( lrot.c, lrot.a );
}
*/

Projectile3D * Turret::fireProjectile3D( const Vec3d& pos0, const Mat3d& rot0, const Vec3d& gvel  ){
    //printf( " Gun fireProjectile3D \n" );
	Projectile3D * p = new Projectile3D();
	//Mat3d rotmat;
	//globalPos( pos0, rot0, p->pos  );
	//globalRot( rot0,        rotmat );
	p->pos.set( gpos );
	p->vel.set_mul( grot.a, muzzle_velocity );
	p->vel.add( gvel );
	p->setMass( Projectile3D_mass );
	//p->update_old_pos();
    //printf( " Gun fireProjectile3D DONE \n" );
	return p;
}
