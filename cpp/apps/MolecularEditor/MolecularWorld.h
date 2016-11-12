
#ifndef MolecularWorld_h
#define MolecularWorld_h

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"

#include "AtomTypes.h"
#include "MoleculeType.h"

class MolecularWorld{
	public:
	int nmols;
	int nAtomTypes,nMolTypes;
	MoleculeType *  molTypes;
	MoleculeType ** instances;
	Vec3d  *pos,*vpos,*fpos;
	Quat4d *rot,*vrot,*frot;

	// temporary variables
	int nptmp;
	Vec3d  *Tps_i, *Tps_j;
	Vec3d  *fs_i,  *fs_j;
	int nInteractions;

	AtomTypes atomTypes;
	double    *C6s,*C12s;

	double surf_z0;
	double surf_zMin;
	double surf_Emin;
	Vec3d  surf_hat;

	DynamicOpt * optimizer;

// ======== initialization

	void initParams( );
	void initTPoints();
	int  loadMolTypes ( char const* dirName, char const* fileName );
	int  loadInstances( char const* filename     );
	bool fromDir      ( char const* dirName, char const* atom_fname, char const* mol_fname, char const* instance_fname );

    inline MolecularWorld(){};
	//MolecularWorld( char const* filename, MoleculeType * molTypeList );
	//MolecularWorld( int nmols_, MoleculeType ** molecules_ );
	void makeFF( );

// =========== Rotation optimization

	void transformPoints( const Vec3d& pos, const Quat4d& rot, int npoints, Vec3d * points, Vec3d * Tpoints );
	void forceFromPoints( int npoints, Vec3d * points, Vec3d * forces,  const Quat4d& q,  Vec3d& fp, Quat4d& fq );

	void cleanPointForce( int npoints, Vec3d * forces );
	void assembleForces( );

	void rigidOptStep( );

// =========== view utils

	//void draw(){

};

#endif
