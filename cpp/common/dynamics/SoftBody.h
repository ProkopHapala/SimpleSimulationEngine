
#ifndef SoftBody_h
#define SoftBody_h

#include "fastmath.h"
#include "Vec3.h"

const int def_nMaterials = 4;                            //  Steel    Dural   Ti     CFRP
static double def_materials_density [def_nMaterials] = { 7.8e+3, 2.7e+3, 4.4e+3, 1.7e+3 };
//static double def_bondTypes_kTens   [def_nMaterials] = { 200e+9,  69e+9, 110e+9,  50e+9 };
//static double def_bondTypes_kPres   [def_nMaterials] = { 200e+9,  69e+9, 110e+9,  50e+9 };

static double def_bondTypes_kTens   [def_nMaterials] = { 200e+8,  69e+8, 110e+8,  50e+8 };
static double def_bondTypes_kPres   [def_nMaterials] = { 200e+8,  69e+8, 110e+8,  50e+8 };

const int     def_nBondTypes = 3;                         // steel, graphene,
static double def_bondTypes_area     [def_nBondTypes]={ 1e-4,   1e-4,   1e-4 }; // [ m^2 ]
static int    def_bondTypes_material [def_nBondTypes]={    1,      1,      1 };
/*
static double def_bondTypes_density[def_nBondTypes]={ 7.8e+3, 7.8e+3, 7.8e+3 }; // [ kg/m^3 ]
static double def_bondTypes_kTens  [def_nBondTypes]={ 200e+9, 200e+9, 200e+9 }; // [ Pa ]
static double def_bondTypes_kPres [def_nBondTypes]={ 200e+9, 200e+9, 200e+9 }; // [ Pa ]
*/

static int def_fix[6] = { 0, 1, 2, 3, 4, 5 };

// ==================
//    BondTypes
// ==================

class BondTypes{
	public:
	int n;
	double * area;
	double * density;
	double * kTens;
	double * kPres;

	BondTypes( int n_, double * area_, double * density_, double * kTens_,  double * kPres_ ){
		n = n_;
		area    = area_; density = density_; kTens   = kTens_; kPres  = kPres_;
	};

	BondTypes( int n_, double * area_, int * materialBook ){
		n = n_;
		area    = area_;
		density = new double[n];
		kTens   = new double[n];
		kPres  = new double[n];
		for (int i=0; i<n; i++){
			int imat = materialBook[i];
			density[i] = def_materials_density[imat];
			kTens  [i] = def_bondTypes_kTens  [imat];
			kPres  [i] = def_bondTypes_kPres  [imat];
		}
	};
};

// ==================
//    SoftBody
// ==================

class SoftBody{
	public:
	// points
	int npoints;
	double * mass;
	double * drag;
	Vec3d  * points;
	// axuliary
	Vec3d  * velocities;
	Vec3d  * forces;
	double * invMass;

	// bonds
	int nbonds;
	int    * bonds;
	double * kTens;
	double * kPres;
	double * l0s;

	// constrains
	int nfix;
	int * fix;

	// ==== function declarations


	int insertBond( int i, int j, double kTens_, double kPres_, double l0s_ ){
	    int ib2 = nbonds<<1;
        bonds[ib2  ] = i;
        bonds[ib2+1] = j;
        kTens[i] = kTens_;
        kPres[i] = kPres_;
        l0s  [i] = l0s_;
	};




	void bondFromType( int ib, int * bondTypes, const BondTypes& bondTypesBooks );
	void evalForces( const Vec3d& gravity, const Vec3d& airFlow );
	void applyConstrains();
	void move( double dt, double damp );
	void draw( float forceScale );


	void init(
		int npoints_, int nbonds_, int nfix_,
		Vec3d  * points_, double * mass_, double * drag_,
		int    * bonds_, int * bondTypes, const BondTypes& bondTypeBooks,
		int    * fix_
	);

	void init(
		int npoints_, int nbonds_, int nfix_,
		Vec3d  * points_, double * mass_,    double * drag_,
		int    * bonds_,  double * kTens_,   double * kPres_,   double * l0s_,
		int    * fix_
	);

	// ===== inline functions

	inline double getBondLength( int ib ){
		int ib2 = ib<<1;
		int i   = bonds[ib2  ];
		int j   = bonds[ib2+1];
		Vec3d d; d.set_sub( points[i], points[j] );
		return  d.norm();
	}

	inline void evalBondForce( int ib ){
		int ib2 = ib<<1;
		int i   = bonds[ib2  ];
		int j   = bonds[ib2+1];
		Vec3d d; d.set_sub( points[i], points[j] );
		double l  = d.norm(); // this should be optimized
		double dl = ( l - l0s[ib] ) / l;
		if( dl > 0 ){
			d.mul( kTens[ib]*dl );
		}else{
			d.mul( kPres[ib]*dl );
		}
		forces[j].add( d );
		forces[i].sub( d );
	}

	inline void evalPointForce( int i, const Vec3d& gravity, const Vec3d& airFlow ){
		forces[i].set_mul( gravity, mass[i] ); // here we clear forces
		Vec3d vrel; vrel.set_sub( airFlow, velocities[i] );
		forces[i].add_mul( vrel, drag[i] * vrel.norm() );
	}

};

#endif

