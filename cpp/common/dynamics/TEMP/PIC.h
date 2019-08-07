
#ifndef PIC_h
#define PIC_h

#include "fastmath.h"
#include "Vec3.h"
#include "Grid3D.h"
//#include "DynamicOpt.h"

/*
 Particle in cell simulation of visco-elastic material with intendent application in:
  - simulation of projectile impact on armor
  - shaped charges
  - explosions and ablation acceleration - including nuclear
*/





class SoftBody{
	public:
	// points
	int      npoints;
	Vec3d  * points     = NULL;
	Vec3d  * velocities = NULL;
	Vec3d  * forces     = NULL;
    // parameters
    double * mass       = NULL;
	double * drag       = NULL;
	double * invMass    = NULL;

	// bonds
	int    nbonds;
    Bond *  bonds       = NULL;

	// constrains
	int   nfix=0;
	int * fix = NULL;

	bool own_points, own_mass, own_fix;

    Grid3D grid;
    
	// ==== function declarations


    void update_neighbors(  );
	void evalForces      (  );
	void applyConstrains (  );
	void move_LeapFrog   (  );
	void step            (  );

    void deallocateAll( );
    void allocate( int npoints_, int nbonds_, int nfix_ );
    void setPoints( int npoints_,  Vec3d  * points_, double * mass_, double * drag_ );
    void setConstrains( int nfix_, int  * fix_  );
    void setBonds     ( int n, int * ips, int * its, BondType * bts );
    int  findBonds( double lmax, BondType * bt );
    void prepareBonds ( bool l0_fromPos );
    void preparePoints( bool clearVelocity, double constDrag, double constMass );



};

#endif

