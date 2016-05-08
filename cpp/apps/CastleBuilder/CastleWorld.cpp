
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include "TerrainHydraulics.h"

#include "CastleWorld.h" // THE HEADER

// ======================  TestApp

class CastleWorld: public TerrainSimplex{
	public:

};

/*
CastleWorld::CastleWorld( ){

    terrain.allocate( 512, 512 );
    //terrain.genTerrainNoise( 14, 0.3, 0.7, 0.6, 45454, {1000.0,1000.0} );
    //terrain.genTerrainNoise( 14, 0.5, 1.0,  0.7, 1.2, 45454, {100.0,100.0} );
    terrain.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );

    //for (int i=0; i<terrain.ntot; i++){  terrain.water[i] = terrain.ground[i]; }

    //for (int i=0; i<terrain.ntot; i++){  terrain.water[i] = 1.0; terrain.water_[i] = 1.0; }
    terrain.initErrosion( 1.0 );

    shape=glGenLists(1);
}

*/













