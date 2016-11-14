
#ifndef MoleculeType_h
#define MoleculeType_h

#include "Vec3.h"

#include "AtomTypes.h"

class MoleculeType{
	public:
	AtomTypes * typeList = NULL;
	int natoms=0,nbonds=0;
	int    * atypes = NULL;
	Vec3d  * xyzs   = NULL;
	double * Qs     = NULL;
	int    * bonds  = NULL;
	int viewlist=0;

	double Rmax = 0.0;

	// ---- Functions

	void allocateAtoms( int n );

	inline MoleculeType(      ){};
	MoleculeType( int natoms_ );

	// --- function implementation

	bool loadFromFile_bas( char const* filename );
	MoleculeType( char const* filename );

	int findBonds( double sc );

	//int drawAtom( int i, int nsphere, float atomscale, Uint32 color );
	//int drawBond( int i, int j, int nstick, float bondwidth  );
	//int makeViewCPK ( int nsphere, int nstick, float atomscale, float bondwidth );

	void initBounds();
	void toCOG_minmax();
	void toCOG_average();

};

#endif
