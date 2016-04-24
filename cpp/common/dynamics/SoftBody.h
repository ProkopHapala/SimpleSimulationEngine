
#ifndef SoftBody_h
#define SoftBody_h

#include "fastmath.h"
#include "Vec3.h"
//#include "DynamicOpt.h"

// ==================
//    BondTypes
// ==================

class BondType{
	public:
	int id;
	double linearDensity;
    double kPress,kTens;  // stiffness
	double sPress,sTens;  // strength
};

class Bond{
    public:
    uint16_t  id;        // unique indentifier
    uint16_t  i,j;       // end node index
    //double mass;
    bool   broken;
	double l0;           // relaxed length
	BondType type;

    inline double getMass(){ return l0 * type.linearDensity; }
    inline double getDrag(){ return l0 ; }  // this could be improved later

	inline double evalFoce( double l ){
        double dl = ( l - l0 ) / l;
        double f;
        if( dl > 0 ){
            f = type.kTens *dl;
        }else{
            f = type.kPress*dl;
        }
        //printf( "%f %f %f %f\n",    l, l0, dl, f );
        return f;
	}

    inline double evalFoceBreak( double l ){
        double dl = ( l - l0 ) / l;
        double f;
        if( dl > 0 ){
            f = type.kTens*dl;
            if( f >  type.sTens ){
                broken = true;
                return 0;
            }
        }else{
            f = type.kPress*dl;
            if( f >  type.sPress ){
                broken = true;
                return 0;
            }
        }
        return f;
	}

};

// ==================
//    SoftBody
// ==================

class SoftBody{
	public:
	// points
	int npoints;
	Vec3d  * points     = NULL;
	Vec3d  * velocities = NULL;
	Vec3d  * forces     = NULL;
    // parameters
    double * mass     = NULL;
	double * drag     = NULL;
	double * invMass  = NULL;

	// bonds
	int nbonds;
    Bond * bonds;

	// constrains
	int   nfix;
	int * fix;

	bool own_points, own_mass, own_fix;

	//Vec3d gravity = {0.0,-9.81,0.0}, airFlow={0.0,0.0,0.0};
	Vec3d gravity, airFlow;
    double dt   = 0.01;
    double damp = 0.99;

	// ==== function declarations

	void evalForces     (  );
	void applyConstrains(  );
	void move_LeapFrog  (  );
	void step           (  );

    void allocate       ( int npoints_, int nbonds_, int nfix_, Vec3d  * points_, double * mass_, double * drag_,  int * fix_ );
    void prepareBonds ( bool l0_fromPos );
    void preparePoints( bool clearVelocity, double constDrag, double constMass );

	// ===== inline functions

	inline double getBondLength( uint16_t i, uint16_t j ){
		Vec3d d; d.set_sub( points[i], points[j] );
		return d.norm();
	}

	inline void addBondForce( Bond& bond ){
        Vec3d d; d.set_sub( points[bond.i], points[bond.j] );
        double  l = d.norm();      // this should be optimized
        double  f = bond.evalFoce( l );
        //double  f = evalFoceBreak( l );
        d.mul( f );
        //printf( " bond force %i %i %f %f (%3.3f,%3.3f,%3.3f)\n", bond.i, bond.j, l, f, d.x, d.y, d.z );
        forces[bond.j].add( d );
        forces[bond.i].sub( d );
	}

	inline void evalPointForce( int i, const Vec3d& gravity, const Vec3d& airFlow ){
		forces[i].set_mul( gravity, mass[i] ); // here we clear forces
		Vec3d vrel; vrel.set_sub( airFlow, velocities[i] );
		forces[i].add_mul( vrel, drag[i] * vrel.norm() );
	}

};

#endif

