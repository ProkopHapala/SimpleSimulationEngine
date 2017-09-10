#ifndef Body_h
#define Body_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

// ========================
//   CLASS :   KinematicBody
// ========================

class KinematicBody{ public:
	Vec3d lpos = (Vec3d){0.0,0.0,0.0};
	Mat3d lrot = (Mat3d){ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
	inline void globalPos( const Vec3d& pos0, const Mat3d& rot0, Vec3d& gpos ){ rot0.dot_to( lpos, gpos ); gpos.add( pos0 ); }
	inline void globalRot( const Mat3d& rot0,                    Mat3d& grot ){ grot.set_mmul( lrot, rot0 ); }

    inline void globalPosT( const Vec3d& pos0, const Mat3d& rot0, Vec3d& gpos ){ rot0.dot_to_T   ( lpos, gpos ); gpos.add( pos0 ); }
	inline void globalRotT( const Mat3d& rot0,                    Mat3d& grot ){ grot.set_mmul_NT( lrot, rot0 ); }
	//inline void globalRot( const Mat3d& rot0, Mat3d& grot ){ grot.set_mmul_NT( lrot, rot0 ); }
};

// ========================
//   CLASS :   PointBody
// ========================

class PointBody{ public:
	// parameters
	double	mass    = 1.0;
	// auxiliary parameters
	double	invMass = 1.0;
	// State variables
	Vec3d pos = (Vec3d){0.0,0.0,0.0};
	Vec3d vel = (Vec3d){0.0,0.0,0.0};
	// auxiliary variables
	Vec3d force = (Vec3d){0.0,0.0,0.0};

	// ==== function declarations

	virtual void evalForce();
	virtual void move(double dt);
	//virtual void render();

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

// ========================
//   CLASS :   RigidBody
// ========================

class RigidBody : public PointBody { public:
	// parameters
	Mat3d	Ibody    = (Mat3d){ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
	// auxiliary parameters
	Mat3d	invIbody = (Mat3d){ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
	// State variables
	Quat4d qrot  = (Quat4d){0.0,0.0,0.0,1.0};
	Vec3d      L = (Vec3d){0.0,0.0,0.0};
	// auxiliary variables
	Mat3d rotMat = (Mat3d){ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
	Mat3d invI   = (Mat3d){ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
	Vec3d omega  = (Vec3d){0.0,0.0,0.0};
	Vec3d torq   = (Vec3d){0.0,0.0,0.0};

	int shape = 0; // displayList

	// ==== function declarations

	void from_mass_points( int n, double* amass, Vec3d* apos );
	void init( );
	//void apply_anchor( double k, const Vec3d& lpos, const Vec3d& gpos0 );
	virtual void move( double dt );
	//virtual void render();

	// ==== inline functions

	inline void clean_temp( ){  force.set(0.0); torq.set(0.0);  }

	inline void update_aux( ){
        double qr2 = qrot.norm2();
        if( (qr2 > 1.000001d) || (qr2 < 0.999999d) ){
            qrot.mul(1/sqrt(qr2));
        }
		qrot.toMatrix  ( rotMat );
		Mat3d tmp;
		tmp.set_mmul_NT(  invIbody, rotMat  ); invI.set_mmul( rotMat, tmp );
		//tmp.set_mmul(  invIbody, rotMat  ); invI.set_mmul_TN( rotMat, tmp );
		//tmp.set_mmul(  invIbody, rotMat  ); invI.set_mmul_NT( rotMat, tmp );
	};

    inline void move_RigidBody( double dt ){
        // postion
        vel.add_mul( force, dt*invMass );
        pos.add_mul( vel, dt   );
        // rotation
        L   .add_mul    ( torq, dt  );  // we have to use angular momentum as state variable, omega is not conserved
        invI.dot_to     ( L,   omega );
        //invI.dot_to_T     ( L,   omega );
        //qrot.dRot_exact ( dt,  omega );
        qrot.dRot_taylor2( dt,  omega );
        update_aux(   );
        //printf("force (%3.3f,%3.3f,%3.3f) vel (%3.3f,%3.3f,%3.3f) pos (%3.3f,%3.3f,%3.3f)\n", force.x,force.y,force.z, vel.x,vel.y,vel.z,  pos.x, pos.y, pos.z  );
        //printf("L (%3.3f,%3.3f,%3.3f) omega (%3.3f,%3.3f,%3.3f) qrot (%3.3f,%3.3f,%3.3f,%3.3f)\n", L.x,L.y,L.z, omega.x,omega.y,omega.z,  qrot.x, qrot.y, qrot.z, qrot.w  );
    };

    inline void glob2loc( const Vec3d& gp, Vec3d& lp ) const{
        Vec3d tmp; tmp.set_sub(gp,pos);
        //rotMat.dot_to( tmp, lp );
        rotMat.dot_to_T( tmp, lp );
    };

    inline void loc2glob( const Vec3d& lp, Vec3d& gp ) const {
        //rotMat.dot_to_T( lp, gp );
        rotMat.dot_to( lp, gp );
        gp.add(pos);
    };

    inline void velOfPoint( const Vec3d& lp, Vec3d& gv, Vec3d& gdp ) const {
        rotMat.dot_to( lp, gdp    );
        gv.set_cross ( omega, gdp );
        gv.add(vel);
        //gp.add(pos);
    }

	inline void apply_force( const Vec3d& dforce, const Vec3d& gdpos ){
		torq .add_cross( gdpos, dforce );
		//torq .add_cross( dforce, gdpos );
		force.add( dforce );
	};

	inline void apply_anchor( double k, const Vec3d& lpos, const Vec3d& gpos0 ){
		Vec3d sforce, gdpos;
		//rotMat.dot_to_T(  lpos, gdpos   );
		rotMat.dot_to(  lpos, gdpos   );
		sforce.set   (( gdpos + pos - gpos0 )*(-k) );
		apply_force  (  sforce, gdpos );
		//drawLine( gpos0, gdpos + pos  );
	};

	inline void checkStateNormal(){
		// check if rotation is normalized
		double qr2  = qrot.norm2();
		double dqr2 = qr2-1;
		if( (dqr2*dqr2) > 0.0001 ){ qrot.mul( 1/sqrt(qr2) ); }
	}

	inline void initOne(){
	    setMass( 1.0 );
        Ibody.a.set(1,0,0);
        Ibody.b.set(0,1,0);
        Ibody.c.set(0,0,1);
        Ibody.invert_to( invIbody );
        qrot.setOne();
        qrot.toMatrix   ( rotMat );
	};

    inline void setPose( const Vec3d& pos_, const Vec3d& dir, const Vec3d& up ){
        //w->kind = kind; w->id = warriorCount; warriorCount++;
        initOne();
        pos.set           ( pos_    );
        rotMat.a.set      ( dir     );
        rotMat.b.set      ( up      );
        rotMat.c.set_cross( dir, up );
        qrot.fromMatrix   ( rotMat );
        //printf( "pos (%g,%g,%g) qrot (%g,%g,%g,%g)\n", pos.x, pos.x, pos.x, qrot.x,qrot.y,qrot.z,qrot.w );
	}

    inline void initSpherical( double mass, double I ){
	    setMass( mass );
        Ibody.a.set(I,0,0);
        Ibody.b.set(0,I,0);
        Ibody.c.set(0,0,I);
        Ibody.invert_to( invIbody );
        qrot.setOne();
        qrot.toMatrix   ( rotMat );
	};

};

// ===============================
//   CLASS :   SpringConstrain
// ===============================

class SpringConstrain{ public:
	Vec3d     p1,p2;
	RigidBody *b1=NULL,*b2=NULL;
	double kPull,kPush,L0;

	// ==== functiopn declarations
    inline void getPoints( Vec3d& gp1, Vec3d& gp2 ){
        b1->loc2glob(p1,gp1);
		if(b2){ b2->loc2glob(p2,gp2); }else{gp2.set(p2);}
    }

    inline Vec3d getForce( const Vec3d& gp1, const Vec3d& gp2 ){
        Vec3d dp; dp.set_sub(gp2, gp1);
        double r = dp.norm();
        //printf("(%3.3f,%3.3f,%3.3f) %f \n", dp.x, dp.y, dp.z, r);
        if( r > L0){
            dp.mul( kPull*(r-L0)/(r+1e-16) );
        }else{
            dp.mul( kPush*(L0-r)/(r+1e-16) );
        }
        //printf("(%3.3f,%3.3f,%3.3f) \n", dp.x, dp.y, dp.z);
        return dp;
    }

	inline Vec3d apply(){
		Vec3d gp1,gp2;
        getPoints( gp1, gp2 );
        Vec3d f = getForce( gp1, gp2 );
        b1      ->apply_force( f   , gp1-b1->pos );
        if(b2)b2->apply_force( f*-1, gp2-b2->pos );
        return f;
	};
	//void render();
	//SpringConstrain( double k_, RigidBody* b1_, RigidBody* b2_, const Vec3d& p1_, const Vec3d& p2_ );
	SpringConstrain( double kPull_, double kPush_, double L0_, RigidBody* b1_, RigidBody* b2_, const Vec3d& p1_, const Vec3d& p2_ );

};

#endif  // #ifndef Body_h
