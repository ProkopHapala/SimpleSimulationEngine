
// according to
// An Introduction to Physically Based Modeling:
// Rigid Body Simulation I—Unconstrained Rigid
// Body Dynamics
// David Baraff
// https://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf

//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

//#include "Draw3D.h"

#include "Body.h" // THE HEADER

// ========================
//   CLASS :   PointBody
// ========================

void PointBody::evalForce()    { force.set( 0.0,-9.81f,0.0 ); };
void PointBody::move(double dt){ move_PointBody(dt);          };
//void PointBody::render()       { Draw3D::drawPoint( pos );    };

// ========================
//   CLASS :   RigidBody
// ========================

void RigidBody::move( double dt ){
    move_RigidBody(dt);
};

void RigidBody::from_mass_points( int n, double* amass, Vec3d* apos ){
	//printf( " a1 \n" );
	mass = 0;
	pos.set(0.0);
	//printf( " a1.1 \n" );
	for(int i=0;  i<n; i++){
		//printf( " %f %f %f \n", apos[i].x, apos[i].y, apos[i].z );
		pos  .add_mul( apos[i], amass[i] );
		mass +=                amass[i];
	};
	//printf( " a2 \n" );
	invMass = 1/mass;
	pos.mul( invMass );
	//printf( " a3 \n" );
	Mat3d Ibody = Mat3dZero;
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
	//printf( " a4 \n" );
	Ibody.invert_to( invIbody );
	//printf( " a5 \n" );
};

/*
void RigidBody::render(){
	glPushMatrix();
	float glmat[16];
	Draw3D::toGLMat( pos, rotMat, glmat );
	glMultMatrixf( glmat );
	glCallList( shape );
	glPopMatrix();
};
*/

void RigidBody::init( ){
	clean_temp( );
	//qrot.toMatrix   ( rotMat );
	//Mat3d tmp; tmp.set_mmul_NT(  invIbody, rotMat  ); invI.set_mmul( rotMat, tmp );
}

// ===============================
//   CLASS :   SpringConstrain
// ===============================


/*
void SpringConstrain::render(){
	Vec3d gp1; b1->rotMat.dot_to( p1, gp1 ); gp1.add( b1->pos );
	Vec3d gp2; b2->rotMat.dot_to( p2, gp2 ); gp2.add( b2->pos );
	Draw3D::drawLine( gp1, gp2 );
}
*/

SpringConstrain::SpringConstrain( double kPull_, double kPush_, double L0_, RigidBody* b1_, RigidBody* b2_, const Vec3d& p1_, const Vec3d& p2_ ){
	kPull=kPull_; kPush=kPush_; L0=L0_; b1=b1_; b2=b2_; p1.set( p1_ ); p2.set( p2_ );
}

