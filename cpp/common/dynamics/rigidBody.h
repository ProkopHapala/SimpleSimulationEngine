
// according to 
// An Introduction to Physically Based Modeling:
// Rigid Body Simulation I—Unconstrained Rigid
// Body Dynamics
// David Baraff
// https://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf

///////////////////////////
//   CLASS :   PointBody
///////////////////////////

class PointBody{
	public:
	// parameters
	double	mass;
	// auxiliary parameters
	double	invMass; 
	// State variables 
	Vec3d pos;
	Vec3d vel;
	// auxiliary variables
	Vec3d force;

	inline void move_PointBody( double dt ){
		vel.add_mul( force, dt*invMass );
		pos.add_mul( vel,   dt );
		//printf( "dt: %f force: ",dt ); printVec( force ); printf( " vel: " ); printVec( vel ); printf( " pos: " ); printVec( pos ); printf( "\n" );
	};

/*
	inline void move_PointBody_fma( double dt ){
		vel.add_mul( force, dt );
		pos.add_mul( vel,   dt );
	};
*/

	inline void clean_temp( ){  force.set(0.0); }

	virtual void evalForce()   { force.set( 0.0,-9.81f,0.0 ); }
	virtual void move(double dt){ move_PointBody(dt);      }

	virtual void render(){
		drawPoint( pos );
	}
};

//////////////////////////////
//   CLASS :   KinematicBody
/////////////////////////////

class KinematicBody{
	public:
	Vec3d lpos;
	Mat3d lrot;
	void globalPos( const Vec3d& pos0, Vec3d& gpos ){ gpos.set_add ( lpos, pos0 ); }
	void globalRot( const Mat3d& rot0, Mat3d& grot ){ grot.set_mmul( lrot, rot0 ); }
};


///////////////////////////
//   CLASS :   RigidBody
///////////////////////////


class RigidBody : public PointBody {
	public:
	// parameters
	Mat3d	Ibody;
	// auxiliary parameters
	Mat3d	invIbody; 
	// State variables        
	Quat4d qrot;             
	Vec3d      L;            
	// auxiliary variables
	Mat3d rotMat;           
	Mat3d invI;        
	Vec3d omega; 
	Vec3d torq; 

	int shape; // displayList 

	inline void clean_temp( ){  force.set(0.0); torq.set(0.0);  }

	inline void update_aux( ){
		qrot.toMatrix   ( rotMat );
		Mat3d tmp; tmp.set_mmul_NT(  invIbody, rotMat  ); invI.set_mmul( rotMat, tmp ); 
	};

	inline void move( double dt ){
		// postion
		vel.add_mul( force, dt*invMass );
		pos.add_mul( vel, dt   );
		// rotation
		L   .add_mul    ( torq, dt  );  // we have to use angular momentum as state variable, omega is not conserved
		invI.dot_to     ( L,   omega );
		//qrot.dRot_exact ( dt,  omega );
		qrot.dRot_taylor2( dt,  omega );
		update_aux(   );
	};

	inline void apply_force( const Vec3d& dforce, const Vec3d& gdpos ){
		torq  .add_cross( gdpos, dforce );
		force.add( dforce ); 
	};

	inline void apply_anchor( double k, const Vec3d& lpos, const Vec3d& gpos0 ){
		Vec3d sforce, gdpos;
		rotMat.dot_to(  lpos, gdpos   );
		sforce.set   (( gdpos + pos - gpos0 )*(-k) );
		apply_force  (  sforce, gdpos );
		drawLine( gpos0, gdpos + pos  );
	};

	void from_mass_points( int n, double* amass, Vec3d* apos ){
		printf( " a1 \n" );
		mass = 0;
		pos.set(0.0);
		printf( " a1.1 \n" );
		for(int i=0;  i<n; i++){
			printf( " %f %f %f \n", apos[i].x, apos[i].y, apos[i].z );
			pos .add_mul( apos[i], amass[i] );  
			mass +=                amass[i];
		};
		printf( " a2 \n" );
		invMass = 1/mass;
		pos.mul( invMass );
		printf( " a3 \n" );
		for(int i=0;  i<n; i++){
			double mi = amass[i];
			Vec3d d; d.set( apos[i] - pos );
			double xx = d.x*d.x; double yy = d.y*d.y;  double zz = d.z*d.z; 
      		double xy = d.x*d.y; double xz = d.x*d.z;  double yz = d.y*d.z; 
			Ibody.xx += mi*(  yy + zz );
      		Ibody.xy += mi*( -xy );
      		Ibody.xz += mi*( -xz );
      		Ibody.yz += mi*( -xy );
      		Ibody.yy += mi*(  xx + zz );
      		Ibody.yz += mi*( -yz );
      		Ibody.zx += mi*( -xz );
      		Ibody.zy += mi*( -yz  );
      		Ibody.zz += mi*(  xx + yy  );      
		};
		printf( " a4 \n" );
		Ibody.invert_to( invIbody );
		printf( " a5 \n" );
	};


	void render(){
		glPushMatrix();
		float glmat[16];
		toGLMat( pos, rotMat, glmat );
		glMultMatrixf( glmat );
		glCallList( shape );
		glPopMatrix();
	};


	void init( ){
		clean_temp( );
		qrot.toMatrix   ( rotMat );
		Mat3d tmp; tmp.set_mmul_NT(  invIbody, rotMat  ); invI.set_mmul( rotMat, tmp ); 
	}

};


////////////////////////////////
//   CLASS :   SpringConstrain
////////////////////////////////

class SpringConstrain{
	public:
	Vec3d p1,p2;
	RigidBody *b1,*b2;
	double k;
 
	SpringConstrain( double k_, RigidBody* b1_, RigidBody* b2_, const Vec3d& p1_, const Vec3d& p2_ ){
		k=k_; 
		b1=b1_; 
		b2=b2_;
		p1.set( p1_ );
		p2.set( p2_ );
	}

	void apply(){
		Vec3d gp1; b1->rotMat.dot_to( p1, gp1 );
		Vec3d gp2; b2->rotMat.dot_to( p2, gp2 );
		Vec3d gdp = ( (gp2+b2->pos)-(gp1+b1->pos) );
		b1->apply_force ( gdp*( k), gp1 );
		b2->apply_force ( gdp*(-k), gp2 );
	}

	void render(){
		Vec3d gp1; b1->rotMat.dot_to( p1, gp1 ); gp1.add( b1->pos );
		Vec3d gp2; b2->rotMat.dot_to( p2, gp2 ); gp2.add( b2->pos );
		drawLine( gp1, gp2 );
	}

};

