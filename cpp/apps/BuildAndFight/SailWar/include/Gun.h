class Gun : public KinematicBody {
	public:
	double muzzle_velocity;
	
	Projectile * fireProjectile( const Vec3d& gpos, const Mat3d& grot, const Vec3d& gvel  ){
		Projectile * p = new Projectile();
		Mat3d rotmat;
		globalPos( gpos, p->pos );
		globalRot( grot, rotmat );
		p->vel.set_mul( rotmat.a, muzzle_velocity );  
		p->vel.add( gvel );
	}

};
