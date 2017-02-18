
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

int Turret::shoot( std::vector<Projectile3D*>& projectiles, const Vec3d& gvel ){
    //printf( " Gun fireProjectile3D \n" );
	Projectile3D * p = new Projectile3D();
	//Mat3d rotmat;
	//globalPos( pos0, rot0, p->pos  );
	//globalRot( rot0,        rotmat );
	p->pos.set    ( gpos + gun_dir*type->Rsize );
	p->vel.set_mul( gun_dir, type->vmuzzle );
	p->vel.add    ( gvel );
	p->setMass    ( type->prjmass );
	//p->update_old_pos();
    //printf( " Gun fireProjectile3D DONE \n" );
	//return p;

	printf( "shooting projectile %i (%f,%f,%f) (%f,%f,%f) \n", projectiles.size(), p->pos.x,p->pos.y,p->pos.z,   p->vel.x,p->vel.y,p->vel.z );
	projectiles.push_back(p);

	reload = 0.0;
	return 1;

}

//void Turret::draw(){};
