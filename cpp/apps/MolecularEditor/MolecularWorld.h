
#ifndef MolecularWorld_h
#define MolecularWorld_h

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"

#include "AtomTypes.h"
#include "MoleculeType.h"


class MolecularLink{
    public:
    int i,j;
    Vec3d posi,posj;
    double l0,k;

    inline Vec3d getForce(const Vec3d& dp){
        double r2 = dp.norm2();
        Vec3d f;
        //printf( "r2 %g  l0 %g\n", sqrt(r2), l0 );
        if( r2 > l0*l0 ){
            double r = sqrt(r2);
            f.set_mul( dp, (r-l0)/r );
        }else{
            f.set(0.0);
        }
        return f;
    }
};

class MolecularWorld{
	public:
	int nmols=0,nAtomTypes=0,nMolTypes=0,nLinkers=0;
    AtomTypes atomTypes;
	MoleculeType *   molTypes=NULL;
	MoleculeType ** instances=NULL;
	MolecularLink * linkers=NULL;

	Vec3d  *pos=NULL,*vpos=NULL,*fpos=NULL,*invMpos=NULL;
	Quat4d *rot=NULL,*vrot=NULL,*frot=NULL,*invMrot=NULL;

	bool * constrains = NULL;

	// temporary variables
	int nptmp=0;
	Vec3d  *Tps_i=NULL, *Tps_j=NULL;
	Vec3d  *fs_i =NULL,  *fs_j=NULL;


	double    *C6s=NULL,*C12s=NULL;
	int nInteractions;

	double surf_z0;
	double surf_zMin;
	double surf_Emin;
	Vec3d  surf_hat;

	DynamicOpt * optimizer=NULL;

// ======== initialization

	void initParams( );
	void initTPoints();

	int  loadMolTypes  ( char const* dirName, char const* fileName );
	int  loadInstances ( char const* fileName );
	int  loadLinkers   ( char const* fileName );
	bool fromDir       ( char const* dirName, char const* atom_fname, char const* mol_fname, char const* instance_fname );


	int  exportAtomsXYZ( FILE * pFile, const char * comment );

    inline MolecularWorld(){};
	//MolecularWorld( char const* filename, MoleculeType * molTypeList );
	//MolecularWorld( int nmols_, MoleculeType ** molecules_ );
	void makeFF( );

// =========== Rotation optimization

	void transformPoints( const Vec3d& pos, const Quat4d& rot, int npoints, Vec3d * points, Vec3d * Tpoints );
	void forceFromPoints( int npoints, Vec3d * points, Vec3d * forces,  const Quat4d& q,  Vec3d& fp, Quat4d& fq );

	void cleanPointForce ( int npoints, Vec3d * forces );
	void assembleForces  ( );
	int  applyLinkerForce( );

	void rigidOptStep( );

// =========== view utils

	//void draw(){

};

#endif
