
#ifndef AtomTypes_h
#define AtomTypes_h

#include <stdint.h>

class AtomTypes{
	public:
	int ntypes;
	int      * Zs;
	char    ** names;
	double   * vdwRs;
	double   * vdwEs;
	uint32_t * colors;

	bool loadFromFile( char const* filename );


	inline AtomTypes( ){};
	AtomTypes( char const* filename );
};

#endif



