
class Yacht2D : public RigidBody2D {
	public:
	AeroSurf2D keel;
	AeroSurf2D rudder;
	AeroSurf2D mast;
	
	inline void applySailForces( const Vec2d& windSpeed, const Vec2d& watterSpeed ){
		Vec2d lvel;
		lvel.set_sub( watterSpeed, vel );
		//printf( " keel   "); 
		keel  .assertAeroForce( *this, lvel, 1000.0 );
		//printf( " rudder "); 
		rudder.assertAeroForce( *this, lvel, 1000.0 );
		lvel.set_sub( windSpeed, vel   );
		//printf( " mast   "); 
		mast.assertAeroForce  ( *this, lvel, 1.225  );
		glColor3f( 0.5f, 0.1f, 0.1f ); drawVecInPos( force*0.01, pos );
	}

	virtual void draw( ){ 
		keel.draw  ( *this );
		rudder.draw( *this );
		mast.draw  ( *this );
	}
};





