
class Projectile : public PointBody {
	public:
	GameWorld world;
	double dragCoef;
	double penetration;
	double damage;

	void ground_hit( ){
		printf( " Ground hit %f %f \n", pos.x, pos.y );
	}

	double check_hit(  ){
		if( pos.z < world.ground_level ){ ground_hit( ); }
	}

	void addDragForce( const Vec3d& vwind, Vec3d& aeroForce ){
		Vec3d vair;
		vair.set_sub( vel, world.wind_speed );
		double vr2   = vair.norm2();
		double vr    = sqrt( vr2 );
		aeroForce.add_mul( vair, dragCoef * vr );		
	}

	virtual void evalForce( ){ 
		force.set( 0.0,-9.81f,0.0 );
		addDragForce( world.wind_speed, force ); 
	}
  
};
