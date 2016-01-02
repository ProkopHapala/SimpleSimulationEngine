#ifndef Body_h
#define Body_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

class KinematicBody{
	public:
	Vec3d lpos;
	Mat3d lrot;
	inline void globalPos( const Vec3d& pos0, const Mat3d& rot0, Vec3d& gpos ){ rot0.dot_to( lpos, gpos ); gpos.add( pos0 ); }
	inline void globalRot( const Mat3d& rot0,                    Mat3d& grot ){ grot.set_mmul( lrot, rot0 ); }
	//inline void globalRot( const Mat3d& rot0, Mat3d& grot ){ grot.set_mmul_NT( lrot, rot0 ); }
};

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

	// ==== function declarations

	virtual void evalForce();
	virtual void move(double dt);
	virtual void render();

	// ==== inline functions

	inline void setMass( double mass_ ){
		mass    = mass_;
		invMass = 1 / mass;
	};

	inline void move_PointBody( double dt ){
		vel.add_mul( force, dt*invMass );
		pos.add_mul( vel,   dt );
		//printf( "dt: %f force: ",dt ); printVec( force ); printf( " vel: " ); printVec( vel ); printf( " pos: " ); printVec( pos ); printf( "\n" );
	};

	inline void clean_temp( ){  force.set(0.0); }

};

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

	// ==== function declarations

	void from_mass_points( int n, double* amass, Vec3d* apos );
	void init( );
	//void apply_anchor( double k, const Vec3d& lpos, const Vec3d& gpos0 );
	virtual void move( double dt );
	virtual void render();

	// ==== inline functions

	inline void clean_temp( ){  force.set(0.0); torq.set(0.0);  }

	inline void update_aux( ){
		qrot.toMatrix   ( rotMat );
		Mat3d tmp; tmp.set_mmul_NT(  invIbody, rotMat  ); invI.set_mmul( rotMat, tmp ); 
	};

	inline void apply_force( const Vec3d& dforce, const Vec3d& gdpos ){
		torq .add_cross( gdpos, dforce );
		force.add( dforce ); 
	};

	inline void apply_anchor( double k, const Vec3d& lpos, const Vec3d& gpos0 ){
		Vec3d sforce, gdpos;
		rotMat.dot_to(  lpos, gdpos   );
		sforce.set   (( gdpos + pos - gpos0 )*(-k) );
		apply_force  (  sforce, gdpos );
		//drawLine( gpos0, gdpos + pos  );
	};

};

class SpringConstrain{
	public:
	Vec3d p1,p2;
	RigidBody *b1,*b2;
	double k;
 
	// ==== function declarations

	void apply();
	void render();
	SpringConstrain( double k_, RigidBody* b1_, RigidBody* b2_, const Vec3d& p1_, const Vec3d& p2_ );

};

#endif  // #ifndef Body_h
