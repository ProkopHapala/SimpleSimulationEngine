
///////////////////////////
//   CLASS :   PointBody2D
///////////////////////////

class PointBody2D{
	public:
	// parameters
	double mass;
	// auxiliary parameters
	double invMass;
	// State variables 
	Vec2d pos;
	Vec2d vel;
	// auxiliary variables
	Vec2d force;

	inline void move_PointBody2D( double dt ){
		vel.add_mul( force, dt*invMass );
		pos.add_mul( vel,   dt );
	};

	inline void clean_temp( ){  force.set(0.0); }

	virtual void evalForce()   {  }
	virtual void move( double dt ){ move_PointBody2D(dt);  }
	virtual void draw( ){ drawPointCross( pos, 0.1 );  }

};


///////////////////////////
//   CLASS :   RigidBody2D
///////////////////////////


const double ROT_NORM2_PREC = 1e-6;


class RigidBody2D : public PointBody2D {
	public:
	// parameters
	double	I;
	// auxiliary parameters
	double	invI; 
	// State variables         
	double phi;          
	double omega;            
	// auxiliary variables
	Vec2d   rot;  
	double torq;

	int shape; // displayList 

	inline void clean_temp (             ){ force.set(0.0);          torq=0;                    }
	inline void setDefaults(             ){	omega = 0.0;		     clean_temp();	            }
	inline void setAngle   ( double phi_ ){	phi = phi_;		         rot.fromAngle( phi_ ); 	}
	inline void sinc_rot   (             ){	rot.x = cos( phi );		 rot.y = sin( phi );	}
	inline bool check_rot  (             ){	double r2 = rot.norm2(); return ( fabs( r2 - 1.0d ) > ROT_NORM2_PREC );	}

	inline void move_RigidBody2D( double dt ){
		// postion
		vel.add_mul( force, dt*invMass );
		pos.add_mul( vel,   dt );
		// rotation
		omega       += torq  * invI * dt;
		double dphi  = omega * dt; 
		phi         += dphi;
		//printf( " invI, dt, torq, omega, dphi, phi %f %f %f %f %f \n",  invI, dt, torq, omega, dphi, phi  );
		rot.rotate_taylor2( dphi );
		//rot.rotate_taylor2( dphi );
		//if ( check_rot() ) { sinc_rot(); }
		rot.x = cos(phi);
		rot.y = sin(phi);

		//rot.rotate( dphi );
		//double xphi = cos( phi );
		//double yphi = sin( phi );
		//printf( " dphi, phi, rot.x, rot.y, xphi, yphi  %f %f   %f %f   %f %f \n", dphi, phi, rot.x, rot.y, xphi, yphi );

		printf( " torq, omega, phi %f %f %f  vel %f %f \n", torq, omega, phi, vel.x, vel.y  );
	};

	void from_mass_points( int n, double * amass, Vec2d * apos ){
		mass = 0;
		pos.set(0.0);
		for(int i=0; i<n; i++){
			// printf( " %f %f %f \n", apos[i].x, apos[i].y, apos[i].z );
			pos .add_mul( apos[i], amass[i] );  
			mass +=                amass[i];
		};
		invMass = 1/mass;
		pos.mul( invMass );
		for(int i=0; i<n; i++){
			Vec2d d; d.set( apos[i] - pos ); 
			I += amass[i] * d.norm2();
		};
		invI = 1/I;
	};

	inline void apply_force( const Vec2d& dforce, const Vec2d& gdpos ){
		torq += gdpos.cross( dforce );
		force.add( dforce ); 
	};

	inline void apply_anchor( double k, const Vec2d& lpos, const Vec2d& gpos0 ){
		Vec2d sforce, gdpos;
		gdpos.set_mul_cmplx( rot, lpos );
		sforce.set_add( gdpos, pos );
		sforce.sub( gpos0 );
		sforce.mul( -k );
		apply_force  (  sforce, gdpos );
		//printf( " gdpos %f %f \n", gdpos.x, gdpos.y );
		//drawLine( gpos0, gdpos + pos  );
	};

	virtual void move( double dt ){ move_RigidBody2D(dt);                             }
	virtual void draw(           ){ drawPointCross( pos, 0.1 ); drawVecInPos( rot, pos );  }

	virtual void draw_shape( ){ 
		glPushMatrix();
		//glTranslatef( pos.x, pos.y , 0 );
		//glRotatef( phi*(180/M_PI), 0, 0, 1 );
		float mat[16];
		mat[0 ] = +rot.x; mat[1 ] = +rot.y; mat[2 ] = 0;  mat[3 ] = 0;
		mat[4 ] = -rot.y; mat[5 ] = +rot.x; mat[6 ] = 0;  mat[7 ] = 0;
		mat[8 ] = 0;      mat[9 ] = 0;      mat[10] = 1;  mat[11] = 0;
		mat[12] = pos.x;  mat[13] = pos.y;  mat[14] = 0;  mat[15] = 1;
		glMultMatrixf( mat );
		glCallList( shape ); 
		glPopMatrix();
	}

};


////////////////////////////////
//   CLASS :   SpringConstrain2D
////////////////////////////////

class SpringConstrain2D{
	public:
	Vec2d p1,p2;
	RigidBody2D *b1,*b2;
	double k;
 
	SpringConstrain2D( double k_, RigidBody2D* b1_, RigidBody2D* b2_, const Vec2d& p1_, const Vec2d& p2_ ){
		k=k_; 
		b1=b1_; 
		b2=b2_;
		p1.set( p1_ );
		p2.set( p2_ );
	}

	void apply(){
		Vec2d gp1; gp1.set_mul_cmplx( b1->rot, p1 );
		Vec2d gp2; gp2.set_mul_cmplx( b2->rot, p2 );
		Vec2d gdp; gdp.set(  (gp2+b2->pos)-(gp1+b1->pos) );
		b1->apply_force ( gdp*( k), gp1 );
		b2->apply_force ( gdp*(-k), gp2 );
	}

	void draw(){
		Vec2d gp1; gp1.set_mul_cmplx( p1, b1->rot ); gp1.add( b1->pos );
		Vec2d gp2; gp2.set_mul_cmplx( p2, b2->rot ); gp2.add( b2->pos );
		drawLine( gp1, gp2 );
		drawPointCross( gp1, 0.1 );
		drawPointCross( gp2, 0.1 );
	}

};







