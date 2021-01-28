
#ifndef SoftBody_h
#define SoftBody_h

#include "fastmath.h"
#include "Vec3.h"

#include "Lingebra.h"
//#include "LinearElasticity.h"

//#include "DynamicOpt.h"

// ==================
//    BondTypes
// ==================

class BondType{ public:
	int id;
	double linearDensity;
    double kPress,kTens;  // stiffness
	double sPress,sTens;  // strength

	//BondType()=default;
	//BondType()

	void fromThickness( double thick, double modul, double strength, double dens ){
        double area   = thick*thick;
        double k      = modul    * area;
        double s      = strength * area;
        linearDensity = dens     * area;
        printf( "BondType area %g[m^2] stiffness %g[N/m] strength %g[N] dens %g[kg/m]\n", area, k, s, linearDensity );
        kPress=k; kTens=k;
        sPress=s; sTens=s;
	}

	static BondType stick(int id, double thick, double modul, double strength, double dens ){
        BondType bt; bt.id=id; bt.fromThickness( thick, modul, strength, dens );
        return bt;
	}
	static BondType rope(int id, double thick, double modul, double strength, double dens ){
        BondType bt; bt.id=id; bt.fromThickness( thick, modul, strength, dens );
        bt.kPress=0; bt.sPress=0;
        return bt;
	}

};


class Bond{ public:

    uint16_t  id;        // unique indentifier
    uint16_t  i,j;       // end node index
    //double mass;
    bool    broken;
	double  l0;           // relaxed length
	BondType * type;

    inline double getMass()const{ return l0 * type->linearDensity; }
    inline double getDrag()const{ return l0 ; }  // this could be improved later

	inline double evalFoce( double l )const{
        //double dl = ( l - l0 ) / l;
        double dl = ( l - l0 ) / l0;
        //double inv_l0 = 1./l0;
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
        //double dl = ( l - l0 ) / l;
        double dl = ( l - l0 ) / l0;
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

    inline double enforce( double l, double fmax )const{
        double dl     = ( l - l0 ) / l;
        double f;
        if( dl > 0 ){ f = type->kTens  * dl;   }
        else        { f = type->kPress *-dl;  }
        if( f>fmax){ return dl; }else{ return 0; };
        return 0;
	}

};


/*
class BondLinearized{
    uint16_t  i,j;
    Vec3d   dir;
    double  f0;   // force at reference position (pre-strain)
    double  kPress,kPull;

    inline void linearize( Vec3d* ps, double l0 ){
        Vec3d  d = ps[j] - ps[j];
        double l = d.normalize();
        dir[i] = d;

        f0[i]  = (l-l0)*k0;
    }

    inline void addForce( Vec3d* disps, Vec3d* forces ) const {
            Vec3d hat  = dirs[il];
            double dfl = hat.dot( disps[i] - disps[j] );
            if( dfl>0 ){ dfl*=kPress; }else{ dfl*=kPull; }
            hat.mul(dfl);  // f = k * <h|di-dj> * h
            fs[ij.a].add(hat);
            fs[ij.b].sub(hat);
    }

}
*/


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

class SoftBody{ public:
	// points
	int npoints;
	Vec3d  * points     = NULL;
	Vec3d  * velocities = NULL;
	Vec3d  * forces     = NULL;
    // parameters
    double * mass       = NULL;
	double * drag       = NULL;
	double * invMass    = NULL;

	// linearized Sticks
	Vec3d  * disps = NULL;
	Vec3d  * dirs  = NULL;
    double * f0s   = NULL;

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

	//void evalForces     (  );
	void cleanForces    (  );
	void evalBondForces (  );
	void evalPointForces(  );
	void applyConstrains(  );
	void move_LeapFrog  (  );
	void relaxStepGS( double fmax );
	void step           (  );


	// linearized Sticks
	void evalForceLinearizedBonds( );
	void linearizedBonds( );
	void disp2pos( );

    void deallocateAll( );
    void allocate( int npoints_, int nbonds_, int nfix_=0 );
    void setPoints( int npoints_,  Vec3d  * points_, double * mass_, double * drag_ );
    void setConstrains( int nfix_, int  * fix_  );
    void setBonds     ( int n, int * ips, int * its, BondType * bts );
    int  findBonds( double lmax, BondType * bt );
    void prepareBonds ( bool l0_fromPos );
    void preparePoints( bool clearVelocity, double constDrag, double constMass );

	// ===== inline functions

	inline double getBondLength( uint16_t i, uint16_t j )const{
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

    inline void addBondForceLinear( Bond& bond ){
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
		forces[i].add_mul( gravity, mass[i] ); // here we clear forces
		if( viscosity > 0.0 ){
            Vec3d vrel; vrel.set_sub( airFlow, velocities[i] );
            forces[i].add_mul( vrel, viscosity * drag[i] * vrel.norm() );
		}
	}

};


// ==================
//    SoftBodyLinearized
// ==================

class SoftBodyLinearized : public LinSolver { public:

    int  npoints;
    Vec3d  * poss    = NULL;
    Vec3d  * Fextern = NULL;
    //Vec3d  * Fwork   = NULL;
    double * anchorKs = NULL;

    int  nsticks;
    Vec3d  * disps = NULL;  // displacements
    Vec2i  * ijs   = NULL;  // links
    double * l0s   = NULL;  // stick neutral length
    double * ks    = NULL;  // stick stiffness
    Vec3d  * dirs  = NULL;  // stick normalized direction
    //Vec3d  * kDirs;

    void init(int nsticks_, int npoints_, Vec3d* poss_, Vec2i* ijs_, double* ks_=0 ){
        //_realloc(npoints);
        //_realloc(ijs,nsticks);
        nsticks=nsticks_;
        npoints=npoints_;
        poss=poss_;
        ijs=ijs_;

        //LinSolver::init( npoints*3 );
        _realloc(dirs,    nsticks);
        _realloc(disps,   npoints);
        //_realloc(Fwork,   npoints);
        _realloc(Fextern, npoints);
        _realloc(anchorKs, npoints);  for(int i=0; i<npoints; i++){ anchorKs[i]=0.0; }

        if(ks_==0){
            _realloc(ks,nsticks);
            for(int i=0; i<nsticks; i++){ ks[i]=1.0; }
        }else{ ks=ks_; }
        //l0s = l0s_;
    }

    void prepareSticks( bool bSetL0s ){
        for( int il=0; il<nsticks; il++  ){
            const Vec2i& ij    = ijs[il];
            Vec3d d     = poss[ij.a] - poss[ij.b];
            double l    = d.norm();
            d.mul( 1/l ); // displacement_ij = ( pos_i - pos_j )/f
            dirs[il] = d;
            //printf( " %i -> %i,%i  %f,%f,%f   \n", il, ij.a, ij.b, d.x, d.y, d.z );
            if(l0s){
                if(bSetL0s){
                    l0s[il] = l;
                }else{
                    double f0   = ks[il]*(l-l0s[il]);
                    d.mul( f0 );
                    Fextern[ij.a].add(d); // forces due to pre-strain; external force to keep stick under given strain
                    Fextern[ij.b].sub(d);
                }
            }
        }
    }

    void prepare( bool bSetL0s ){
        // - evaluate stick lengths and normalized directions
        for(int i=0; i<npoints; i++ ){
            disps[i]  .set(0.0);
            Fextern[i].set(0.0);
        }
        prepareSticks(true);
        setLinearProblem( npoints*3, (double*)disps, (double*)Fextern, 0 );
    }

    void move(double h){
        for(int i=0; i<npoints; i++){
            poss[i].add_mul( disps[i], h );
            disps[i].set(0.0);
        };
        prepareSticks(false);
    }

    //void disp2force( int nds, int nfs, double * ds_, double * fs_ ){
    void disp2force( int n, Vec3d* ds, Vec3d* fs ){
        //int n=nds/3;
        //Vec3d * ds = (Vec3d*)ds;
        //Vec3d * fs = (Vec3d*)fs;
        //for( int i=0; i<nfs; i++ ){ fs_[i]=0; }
        //for( int i=0; i<n; i++ ){ fs[i].set(0.0); }
        //printf("DEBUG 1.0 \n");
        //for( int i=0; i<npoints; i++ ){ fs[i] = Fextern[i]; }
        for( int i=0; i<npoints; i++ ){
            //fs[i].set(0.0);
            fs[i].set_mul( ds[i], anchorKs[i] );
        }
        //printf("DEBUG 1.1 \n");
        for( int il=0; il<nsticks; il++ ){
            Vec2i ij   = ijs [il];
            Vec3d hat  = dirs[il];
            double dfl = ks  [il] * hat.dot( ds[ij.a] - ds[ij.b] );
            hat.mul(dfl);  // f = k * <h|di-dj> * h
            //Vec3d f    = dirs[il] * dfl;
            //Vec3d f   = kDirs[il] * ( ds[ij.a] - ds[ij.b] );
            fs[ij.a].add(hat);
            fs[ij.b].sub(hat);
            //printf( " %i   %i,%i    %f,%f    %f,%f,%f   %f,%f,%f   %f,%f,%f \n", il, ij.a, ij.b,   dfl,ks[il],  hat.x,hat.y,hat.z,    ds[ij.a].x, ds[ij.a].y, ds[ij.a].z,      ds[ij.b].x, ds[ij.b].y, ds[ij.b].z );
            //printf( " %i   %i,%i    %f,%f    %f,%f,%f   %f,%f,%f   %f,%f,%f \n", il, ij.a, ij.b,   dfl,ks[il],  hat.x,hat.y,hat.z,    fs[ij.a].x, fs[ij.a].y, fs[ij.a].z,      fs[ij.b].x, fs[ij.b].y, fs[ij.b].z );
        }
        //printf("DEBUG 1.2 \n");
    };

    virtual void dotFunc( int n, double * x, double * Ax ) override {
        disp2force( n, (Vec3d*)x, (Vec3d*)Ax );
    }

    /*
    void solve(){
        //prepare();
        //Lingebra::genLinSolve_CG<disp2force>( npoints*3, (double*)disps, (double*)Fextern );
        //SoftBodyLinearized* T;
        Lingebra::genLinSolve_CG( npoints*3, (double*)disps, (double*)Fwork,
            [&](int nds, int nfs, double * ds_, double * fs_){ this->disp2force( nds/3, (Vec3d*)ds_, (Vec3d*)fs_ ); }
            //[&](){ this->disp2force( , (Vec3d*)disps, (Vec3d*)Fwork ); }
        );
    }
    */

};



#endif

