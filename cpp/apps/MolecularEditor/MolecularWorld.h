
#ifndef MolecularWorld_h
#define MolecularWorld_h

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"
#include "TileBuffer3D.h"

#include "radial_splines.h"
#include "AtomTypes.h"
#include "MoleculeType.h"

#include "PointCloudComparator.h"

class MolecularLink{
    public:
    int i,j;
    Vec3d posi,posj;
    double l0,k;
};

class MolecularBond{
    public:
    int imol,jmol;
    int iatom,jatom;
    double l0,k,dlmax;
};

class MolecularWorld{
	public:
	int nmols=0,nMolTypes=0,nLinkers=0,nBonds=0,nSplines=0;
	//int nAtomTypes=0;
    AtomTypes atomTypes;
	MoleculeType   * molTypes=NULL;
	MoleculeType  ** instances=NULL;
	MolecularLink  * linkers=NULL;
	MolecularBond  * bonds=NULL;
	GeneralSpline  * splines;

	Vec3d  *pos=NULL,*vpos=NULL,*fpos=NULL,*invMpos=NULL;
	Quat4d *rot=NULL,*vrot=NULL,*frot=NULL,*invMrot=NULL;

	bool * constrains = NULL;

	// temporary variables
	int nptmp=0;
	Vec3d  *Tps_i=NULL, *Tps_j=NULL;
	Vec3d  *fs_i =NULL,  *fs_j=NULL;

	//double    *C6s=NULL,*C12s=NULL;

	double surf_z0;
	double surf_zMin;
	double surf_Emin;
	Vec3d  surf_hat;

	DynamicOpt * optimizer=NULL;
    int    nInteractions,nAtomTot;
	bool   nonCovalent = true;
	double fmax;

    double Rcut = 6.0;
	TileBuffer3D<uint16_t,16,16,16,512> boxbuf;

// ======== initialization

	void initParams ( );
	void initTPoints( );
	void setCutoff  ( double Rcut_);

	int  loadMolTypes  ( char const* dirName, char const* fileName );
	int  loadInstances ( char const* fileName );
	int  loadLinkers   ( char const* fileName );
	int  loadBonds     ( char const* fileName );
	int  loadSplines   ( char const* fileName );
	bool fromDir       ( char const* dirName, char const* atom_fname, char const* mol_fname, char const* instance_fname );

    int  saveInstances ( char const* fileName );
	int  exportAtomsXYZ( FILE * pFile, const char * comment );

    inline MolecularWorld(){};
	//MolecularWorld( char const* filename, MoleculeType * molTypeList );
	//MolecularWorld( int nmols_, MoleculeType ** molecules_ );
	void makeFF( );

// =========== Rotation optimization

	void transformPoints ( const Vec3d& pos, const Quat4d& rot, int npoints, Vec3d * points, Vec3d * Tpoints );
	void forceFromPoints ( int npoints, Vec3d * points, Vec3d * forces,  const Quat4d& q,  Vec3d& fp, Quat4d& fq );

	void cleanPointForce ( int npoints, Vec3d * forces );
	void assembleForces  ( );
	int  nonBondingFroces_N2  ();
	int  nonBondingFroces_bbox();
	int  nonBondingFroces_buf ();
	int  applyLinkerForce( );
	int  applyBondForce  ( );

	int checkBonds( double flmin, double flmax );
	//void initBufBox( double Rmax );

	void rigidOptStep( );

// =========== view utils

	//void draw(){

};

#endif
