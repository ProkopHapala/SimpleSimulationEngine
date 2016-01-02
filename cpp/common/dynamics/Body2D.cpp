
//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "drawMath2D.h"

#include "Body2D.h" // THE HEADER

// ===========================
// ====  CLASS : PointBody2D
// ===========================

	void PointBody2D::evalForce()   {  }
	void PointBody2D::move( double dt ){ move_PointBody2D(dt);  }
	void PointBody2D::draw( ){ Draw2D::drawPointCross_d( pos, 0.1 );  }

// ==========================
// ====  CLASS : RigidBody2D
// ==========================

void RigidBody2D::from_mass_points( int n, double * amass, Vec2d * apos ){
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

void RigidBody2D::move( double dt ){
	move_RigidBody2D(dt);
};

void RigidBody2D::draw(           ){
	Draw2D::drawPointCross_d( pos, 0.1 );
	Draw2D::drawVecInPos_d( rot, pos );
};

void RigidBody2D::draw_shape( ){
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
};

// =================================
// ====  CLASS : SpringConstrain2D
// =================================

void SpringConstrain2D::apply(){
	Vec2d gp1; gp1.set_mul_cmplx( b1->rot, p1 );
	Vec2d gp2; gp2.set_mul_cmplx( b2->rot, p2 );
	Vec2d gdp; gdp.set(  (gp2+b2->pos)-(gp1+b1->pos) );
	b1->apply_force ( gdp*( k), gp1 );
	b2->apply_force ( gdp*(-k), gp2 );
};

void SpringConstrain2D::draw(){
	Vec2d gp1; gp1.set_mul_cmplx( p1, b1->rot ); gp1.add( b1->pos );
	Vec2d gp2; gp2.set_mul_cmplx( p2, b2->rot ); gp2.add( b2->pos );
	Draw2D::drawLine_d( gp1, gp2 );
	Draw2D::drawPointCross_d( gp1, 0.1 );
	Draw2D::drawPointCross_d( gp2, 0.1 );
};

SpringConstrain2D::SpringConstrain2D( double k_, RigidBody2D* b1_, RigidBody2D* b2_, const Vec2d& p1_, const Vec2d& p2_ ){
	k=k_;
	b1=b1_;
	b2=b2_;
	p1.set( p1_ );
	p2.set( p2_ );
};

