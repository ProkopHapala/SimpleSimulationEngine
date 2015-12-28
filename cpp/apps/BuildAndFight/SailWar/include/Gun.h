class Gun : public KinematicBody {
	public:
	double muzzle_velocity;
	
	void set_direction( Vec3d dir ){
		Vec3d up; up.set( 0.0d, 0.0d, 1.0d );
		lrot.a.set( dir );
		lrot.c.set_cross( lrot.a, up     );
		lrot.b.set_cross( lrot.c, lrot.a );
	}

	Projectile * fireProjectile( const Vec3d& gpos, const Mat3d& grot, const Vec3d& gvel  ){
		Projectile * p = new Projectile();
		Mat3d rotmat;
		globalPos( gpos, p->pos );
		globalRot( grot, rotmat );
		p->vel.set_mul( rotmat.a, muzzle_velocity );  
		p->vel.add( gvel );
	}

};
