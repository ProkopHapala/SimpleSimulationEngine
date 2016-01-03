
class Projectile : public PointBody {
	public:
	GameWorld world;
	double dragCoef;
	double penetration;
	double damage;

	Vec3d old_pos;

	void ground_hit( ){
		printf( " Ground hit %f %f \n", pos.x, pos.y );
	}

	bool check_hit(  ){
		bool hitted = pos.z < world.ground_level; 
		if( hitted ){ 
			ground_hit( ); 
		};
		return hitted;
	}

	void addDragForce( const Vec3d& vwind, Vec3d& aeroForce ){
		Vec3d vair;
		vair.set_sub( vel, world.wind_speed );
		double vr2   = vair.norm2();
		double vr    = sqrt( vr2 );
		aeroForce.add_mul( vair, dragCoef * vr );		
	}

	virtual void evalForce( ){ 
		force.set( 0.0, 0.0, -9.81f );
		//addDragForce( world.wind_speed, force ); 
	}

	virtual void move( double dt ){ 
		old_pos.set( pos );
		PointBody::move( dt );
	}

	virtual void draw(){
		glBegin(GL_LINES);
			//glVertex3f( (float)( old_pos.x ), (float)( old_pos.y ), 0 );
			glVertex3f( (float)( 0 ), (float)( 0 ), 0 );
			glVertex3f( (float)(     pos.x ), (float)(     pos.y     ), 0 );
		glEnd();
		//printf( " I'm projectile \n");
		printf( " projectile   pos  %10.5f %10.5f %10.5f vel %10.5f %10.5f %10.5f \n",   pos.x, pos.y, pos.z,   vel.x, vel.y, vel.z );
	}
  
};
