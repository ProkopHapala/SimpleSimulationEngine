#ifndef Body2D_h
#define Body2D_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"

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

	// ==== function declarations

	virtual void evalForce();
	virtual void move( double dt );
	virtual void draw( );

	// ==== inline functions

	inline void setMass( double mass_ ){
		mass    = mass_;
		invMass = 1 / mass;
	};

	inline void move_PointBody2D( double dt ){
		vel.add_mul( force, dt*invMass );
		pos.add_mul( vel,   dt );
	};

	inline void clean_temp( ){  force.set(0.0); }

};

class RigidBody2D : public PointBody2D {
	public:
	static constexpr double ROT_NORM2_PREC = 1e-6;
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

	// ==== function declarations

	void from_mass_points( int n, double * amass, Vec2d * apos );
	virtual void move( double dt );
	virtual void draw(           ); 
	virtual void draw_shape(     );

	// ==== inline functions

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
		rot.rotate_taylor2( dphi );
		rot.x = cos(phi);
		rot.y = sin(phi);
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

};

class SpringConstrain2D{
	public:
	Vec2d p1,p2;
	RigidBody2D *b1,*b2;
	double k;

	// ==== function declarations

	void apply();
	void draw();
	SpringConstrain2D( double k_, RigidBody2D* b1_, RigidBody2D* b2_, const Vec2d& p1_, const Vec2d& p2_ );

};

#endif
