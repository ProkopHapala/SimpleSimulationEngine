
#ifndef MolecularWorld_h
#define MolecularWorld_h

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"
#include "TileBuffer3D.h"
#include "Grid3D.h"
#include "CubicRuler.h"

#include "radial_splines.h"
#include "AtomTypes.h"
#include "MoleculeType.h"

#include "PointCloudComparator.h"
#include "SphereTreeND.h"

#include "RBodyConfDyn.h"

class MolecularLink{
    public:
    int    i,j;
    Vec3d  posi,posj;
    double l0,k;
};

class MolecularBond{
    public:
    int    imol,jmol;
    int    iatom,jatom;
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

	//  global atomic positions - will be more efficient for fiture

	int     nAtoms;
	int   * atoms_type = NULL;  // index of atomType
	int   * atoms_mol  = NULL;   // index of molecular instance to which the atom belongs
	Vec3d * atoms_pos  = NULL;


	//double    *C6s=NULL,*C12s=NULL;

	double surf_z0;
	double surf_zMin;
	double surf_Emin;
	Vec3d  surf_hat;

	DynamicOpt * optimizer=NULL;
    int    nInteractions,nAtomTot;
	bool   nonCovalent = true;
	double fmax;

	double Rcollision = 3.0;
    double Rcut = 6.0;
	TileBuffer3D<uint16_t,16,16,16,512> boxbuf;

	//Grid3D<TILE_vector<int>,int,50,50,50> grid;

    CubicRuler ruler;
    std::unordered_multimap<int_fast64_t,int>  atomsMap;


// ======== initialization

	void initParams ( );
	void initTPoints( );
	void initAtoms  ( );
	void setCutoff  ( double Rcut_);

	int  loadMolTypes  ( char const* dirName, char const* fileName );
	int  loadInstances ( char const* fileName );
	int  loadLinkers   ( char const* fileName );
	int  loadBonds     ( char const* fileName );
	int  loadSplines   ( char const* fileName );
	bool fromDir       ( char const* dirName, char const* atom_fname, char const* mol_fname, char const* instance_fname );

    int  getNAtoms();
    int  getAtomTypes  ( int   * buff );
    int  getAtomPos    ( Vec3d * buff );

    uint32_t atom2map ( int i, int ix, int iy, int iz );
    int      atom2map ( int i, double r );
    int      atoms2map( );
    //int      check_collision( int imol, const Vec3d& pos, const Quat4d& qrot );
    int      collisionForce( int imol, const Vec3d& pos, const Quat4d& qrot,  Vec3d& fpos, Quat4d& frot );

	int  exportAtomsXYZ( FILE * pFile, const char * comment );
	int  saveInstances ( char const* fileName );

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
