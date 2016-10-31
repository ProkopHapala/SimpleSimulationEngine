
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
    bool    broken;
	double  l0;           // relaxed length
	BondType * type;

    inline double getMass(){ return l0 * type->linearDensity; }
    inline double getDrag(){ return l0 ; }  // this could be improved later

	inline double evalFoce( double l ){
        double dl = ( l - l0 ) / l;
        double f;
        if( dl > 0 ){
            f = type->kTens *dl;
        }else{
            f = type->kPress*dl;
        }
        //printf( "%f %f %f %f  k %f %f \n",    l, l0, dl, f,  type->kTens, type->kPress );
        return f;
	}

    inline double evalFoceBreak( double l ){
        double dl = ( l - l0 ) / l;
        double f;
        if( dl > 0 ){
            f = type->kTens*dl;
            if( f >  type->sTens ){
                broken = true;
                return 0;
            }
        }else{
            f = type->kPress*dl;
            if( f >  type->sPress ){
                broken = true;
                return 0;
            }
        }
        return f;
	}

};

/*
BondType default_BondType = {
    0,          // id
    7.8,        // density
    1e+5,1e+5,  // stiffness
    1e+9,1e+9   // strength
};
*/

static const BondType default_BondType = {
    0,          // id
    7.8,        // density
    100,100,      // stiffness
    1e+9,1e+9   // strength
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
    double * mass       = NULL;
	double * drag       = NULL;
	double * invMass    = NULL;

	// bonds
	int nbonds;
    Bond * bonds        = NULL;

	// constrains
	int   nfix=0;
	int * fix = NULL;

	bool own_points, own_mass, own_fix;

	//Vec3d gravity = {0.0,-9.81,0.0}, airFlow={0.0,0.0,0.0};
	Vec3d gravity, airFlow;
    double dt   = 0.01d;
    double damp = 0.5d;
    double damp_stick = 0.5d;
    double viscosity = -1.0;
	// ==== function declarations

	void evalForces     (  );
	void applyConstrains(  );
	void move_LeapFrog  (  );
	void step           (  );

    void deallocateAll( );
    void allocate( int npoints_, int nbonds_, int nfix_ );
    void setPoints( int npoints_,  Vec3d  * points_, double * mass_, double * drag_ );
    void setConstrains( int nfix_, int  * fix_  );
    void setBonds     ( int n, int * ips, int * its, BondType * bts );
    int  findBonds( double lmax, BondType * bt );
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
        /*
        if( damp_stick > 0 ){
            Vec3d dv; dv.set_sub( velocities[bond.i], velocities[bond.j] );
            double vf = dv.dot(d);
            if(vf<0) d.mul( damp_stick );
        };
        */
        //printf( " bond force %i %i %f %f (%3.3f,%3.3f,%3.3f)\n", bond.i, bond.j, l, f, d.x, d.y, d.z );
        forces[bond.j].add( d );
        forces[bond.i].sub( d );
	}

	inline void evalPointForce( int i, const Vec3d& gravity, const Vec3d& airFlow ){
		forces[i].set_mul( gravity, mass[i] ); // here we clear forces
		if( viscosity > 0.0 ){
            Vec3d vrel; vrel.set_sub( airFlow, velocities[i] );
            forces[i].add_mul( vrel, viscosity * drag[i] * vrel.norm() );
		}
	}

};

#endif

