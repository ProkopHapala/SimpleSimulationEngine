
class Body2D{
	public:
	// parameters
	double mass;
	// Axiliary parameters
	double imass;
	// stte variables
	Vec2d pos,vel;
	Vec2d dir;
	// temporary variables
	Vec2d force;

	inline void move( double dt ){
		vel.mul(0.8);
		vel.add_mul( force, dt );
		pos.add_mul( vel,   dt );
	}

};

void drawBody( const Body2D& b ){
	glColor3f(0.8,0.8,0.8); draw_Vec2d( b.pos );
	glColor3f(0.8,0.1,0.1); draw_Segment2d( b.pos, b.pos + b.force*100 );
	glColor3f(0.0,0.5,0.1); draw_Segment2d( b.pos, b.pos + b.vel*10  );
	glColor3f(0.8,0.8,0.8);
	//print( b.force );
}

/*

- formation is composed of bodies connected by springs
- force between neighbors is radial spring with equilibrium ( want to be particular dostacne from each other )
- force between oponents is also radial but only repulsive ( keep distance given by reach of weapon )
- force exerted by soldier toward opponent depend on cos(alfa), alfa is angle between direction to oponent and orientation of soldier
- protection depends also on this angle

*/


/*
inline double histeresis_func( double c ){
  return 1.0d/(1.0d + c*c );
}
*/






class UnitType{
	public:
	double r0       = 1;  // spacing in compact formation
	double weapon_r = 3;  // reach of weapon
};



inline double springForceSize( double r2, double r0sq ){
	return ( 1 - sqrt( r0sq/r2 ) ); // exact harmonic
//	return ( r2 - r0sq );
//	return ( r2 - r0sq )/(2*r2);
//	return ( r2 - r0sq )/(4*r2*r2);
}


class Regiment{
	public:
	int n;
	Body2D * bodies;

	UnitType * type;
	double r0sq    = 1  ;
	double kBend   = 0.1;
	double kStretch= 1  ;

	Regiment( int n_, UnitType* type_ ){
		n = n_;
		type = type_;
		bodies = new Body2D[n];
	};

	void cohesionForces(  ){
	  
		// first stick  
		Vec2d posL; posL.set    ( bodies[0].pos       );
		Vec2d pos;  pos .set    ( bodies[1].pos       );
		Vec2d dL;   dL  .set_sub( pos, posL           );   double r2L  = dL.norm2( );
		Vec2d FL;   FL  .set_mul( dL, kStretch*springForceSize( r2L, r0sq ) );	

		print( FL );

		double irL  = 1.0d/sqrt(r2L);
		bodies[0  ].dir  .set( dL.y*irL, -dL.x*irL );
		bodies[0  ].force.set( FL );
	
		Vec2d dR,FR;
		double r2R;

		for(int i=2; i<n; i++){    
			Vec2d posR; posR .set ( bodies[i].pos );
			dR .set_sub ( posR, pos     );              r2R  = dR .norm2( );
			FR .set_mul ( dR, kStretch*springForceSize( r2R, r0sq ) );
			Vec2d dLR;  dLR.set_add ( posR, posL    );   double irLR = 1.0d/dLR.norm();
	
			bodies[i-1].dir.set( dLR.y*irLR, -dLR.x*irLR );
			bodies[i-1].force.x = FR.x - FL.x - ( pos.x-0.5*(posL.x+posR.x) )*kBend;
			bodies[i-1].force.y = FR.y - FL.y - ( pos.y-0.5*(posL.y+posR.y) )*kBend;

			// move to next
			FL  .set( FR   ); 
			posL.set( pos  );
			pos .set( posR );
		}
	
		double irR  = 1.0d/sqrt(r2R);
		bodies[n-1].dir  .set( dR.y*irR, -dR.y*irR );
		bodies[n-1].force.set( -FR.x, -FR.y );

	};

	void move( double dt ){ 
		for(int i=0; i<n; i++){  
			bodies[i].move( dt);
		}
	};

	void draw( ){
		for(int i=0; i<n; i++){  
			draw_circle( bodies[i].pos, type->weapon_r );
			drawBody( bodies[i] );
			if(i>0){ draw_Segment2d( bodies[i-1].pos, bodies[i].pos ); }
		}
	};

};


void colisionForce( Regiment& A, Regiment& B ){
	Body2D * As = A.bodies;
	Body2D * Bs = B.bodies;
	double wr2A = sq( A.type->weapon_r );
	double wr2B = sq( B.type->weapon_r );
	for (int i=0; i<A.n; i++  ){
		for (int j=0; j<B.n; j++ ){
			Vec2d d; d.set_sub( Bs[j].pos, As[i].pos );
			double r2 = d.norm2();
			if( r2 < wr2A ){ 		}
			if( r2 < wr2B ){ 		}
		}
	}
};

















