
#include <SDL2/SDL_opengl.h>

#include "drawMath.h"

#include "testUtils.h"

#include "GameWorld.h" // THE HEADER

static constexpr int npts     = 4;
static constexpr double poss[npts*2] = { -1.0, 0.0,   0.0, -0.1,   0.0, +0.1,   +1.0, 0.0  };
static constexpr double mass[npts  ] = {  10.0, 50.0, 50.0, 10.0  };

void GameWorld::update( ){

	//long long nt1 = getCPUticks();

	double g  = -9.81;
	myCraft->checkStateNormal();
	//for(int iSubStep=0; iSubStep<PHYS_STEPS_PER_FRAME; iSubStep++ ){
	for( int i=0; i<perFrame; i++ ){
		myCraft->clean_temp();
		myCraft->force.set( { 0, g*myCraft->mass, 0 } );
		myCraft->force.add_mul(   myCraft->rotMat.c, 50.0 );
		myCraft->applyAeroForces( {0,0,0} );
		//myCraft->applyAeroForces( {0,10,0} );
		myCraft->move(dt);
	}

	//long long nt2  = getCPUticks();
	//double perStep = double(nt2-nt1)/PHYS_STEPS_PER_FRAME;
	//tickSum+= perStep; stepSum++;
	//printf( " PERFORMANCE: %f ticks/step ( in %i steps ) average: %f \n", perStep, PHYS_STEPS_PER_FRAME, tickSum/stepSum );

};


void GameWorld::makeAeroCraft(){
	const int len = 5;
	double masses[len]  = { 4,1,1,0.5f,0.5f };
	Vec3d  poss[len]    = { {0,-0.05,1}, {-2,0,0}, {2,0,0}, {0,0,-3}, {0,0.2,-3} };
    printf("init AeroCraft\n");

    myCraft = new AeroCraft();
	myCraft->from_mass_points( len, masses, poss );
	myCraft->qrot.setOne();
	//myCraft->qrot.set(0,0,-0.5,1); myCraft->qrot.normalize();

	printVec(myCraft->pos);
	printf("pos\n");

	myCraft->wingLeft .craft = myCraft;
	myCraft->wingRight.craft = myCraft;
	myCraft->elevator .craft = myCraft;
	myCraft->rudder   .craft = myCraft;

	myCraft->wingLeft .lpos.set( poss[1] - myCraft->pos );
	myCraft->wingRight.lpos.set( poss[2] - myCraft->pos );
	myCraft->elevator .lpos.set( poss[3] - myCraft->pos );
	myCraft->rudder   .lpos.set( poss[4] - myCraft->pos );

/*
	myCraft->wingLeft .lpos.set( poss[1] );
	myCraft->wingRight.lpos.set( poss[2] );
	myCraft->elevator .lpos.set( poss[3] );
	myCraft->rudder   .lpos.set( poss[4] );
*/

	myCraft->wingLeft .lrot.set( { 1,0,0, 0,1,0,  0,0,1  } ); myCraft->wingLeft  .lrot.rotate( -0.1, { 0,0,1 } );
	myCraft->wingRight.lrot.set( { 1,0,0, 0,1,0,  0,0,1  } ); myCraft->wingRight .lrot.rotate( +0.1, { 0,0,1 } );
	myCraft->elevator .lrot.set( { 1,0,0, 0,1,0,  0,0,1  } ); myCraft->elevator .lrot.rotate( +0.2, { 1,0,0 } );
	printf("elevator lrot\n");
	printMat( myCraft->elevator.lrot );

	myCraft->rudder   .lrot.set( { 0,1,0, 1,0,0,  0,0,1  } );

	myCraft->wingLeft .C.set( 0.05, wingLiftDefault, 0.05 );
	myCraft->wingRight.C.set( 0.05, wingLiftDefault, 0.05 );
	myCraft->elevator .C.set( 0.05, 1.0, 0.05 ); myCraft->elevator.C.mul( 0.1 );
	myCraft->rudder   .C.set( 0.05, 1.0, 0.05 ); myCraft->rudder  .C.mul( 0.1 );

	myCraft->L.set(0,0,0);
	myCraft->init( );

	myCraft->vel.set(0,0,0);
	myCraft->pos.set(0,500,0);

	printf("Ibody\n");
	printMat(myCraft->Ibody);
	printf("invIbody\n");
	printMat(myCraft->invIbody);

};


int GameWorld::makeBuildingsGrid( int nx, int ny, float sx, float sy, float cx, float cy,  float min_height, float max_height ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
	for (int ix=-nx; ix<nx; ix++){
		float x = ix*sx;
		for (int iy=-ny; iy<ny; iy++){
			float height = randf() * (max_height-min_height) + min_height;
			float y = iy*sy;
			drawBox( x, x + sx*cx, 0, height, y, y + sy*cy, 0.75f, 0.75f, 0.75f );
		}
	}
	glEndList();
	return( ilist );
}

int GameWorld::makeBuildingsClusters( int nclustest, int nmin, int nmax, float minx, float maxx, float miny, float maxy, float min_dist, float max_dist, float min_size, float max_size, float min_height, float max_height ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
	int nboxes = 0;
	for (int icluster=0; icluster<nclustest; icluster++){
		float x0 = randf()*(maxx-minx) + minx;
		float y0 = randf()*(maxy-miny) + miny;
		float nb = round(randf()*(nmax - nmin)) + nmin;
		for (int ib=0; ib<nb; ib++){
			float height = randf() * (max_height-min_height) + min_height;
			float x  = x0 + randf()*(max_dist-min_dist) + min_dist;
			float y  = y0 + randf()*(max_dist-min_dist) + min_dist;
			float dx = randf()*(max_size-min_size) + min_size;
			float dy = randf()*(max_size-min_size) + min_size;
			drawBox( x-dx, x+dx, 0, height, y-dy, y+dy, 0.75f, 0.75f, 0.75f );
			nboxes++;
		};
	};
	printf(" %i buildings \n", nboxes );
	glEndList();
	return( ilist );
}

void GameWorld::makeEnvironment( float sz ){

	//buildings_shape = makeBuildings( 10, 10, 100, 100, 0.5, 0.5, 50, 100 );
	buildings_shape = makeBuildingsClusters( 30, 3, 10,   -sz,         sz,         -sz,          sz,            0, 500,   20, 100,   10, 100 );
	//buildings_shape = makeBuildingsClusters( 30, 3, 10, -VIEW_DEPTH/2, VIEW_DEPTH/2, -VIEW_DEPTH/2, VIEW_DEPTH/2,    0, 500,   20, 100,   10, 100 );
	//buildings_shape= makeBuildingsClusters( 100, 5, 5, -VIEW_DEPTH/2, VIEW_DEPTH/2, -VIEW_DEPTH/2, VIEW_DEPTH/2,    0, 500,   100, 100,   10, 100 );

	double h0    = 1;
	//float tersz = VIEW_DEPTH/2;
	//terrain   = FieldPatch::makeList( 15, { 0.5,   -tersz,-tersz,h0,   tersz,-tersz,h0,  -tersz,tersz,h0,   tersz,tersz,h0   }   );

/*
	Vec3d p1,p2,p3,p4;
	p1.set(-tersz,-tersz,h0);
	p2.set( tersz,-tersz,h0);
	p3.set(-tersz, tersz,h0);
	p4.set( tersz, tersz,h0);
	terrain     = FieldPatch::makeList( 3, Rect( 0.5d,   p1,  p2, p3, p4 )   );
*/
	terrain_shape     = fieldPatch.makeList( 3, { 0.5d,   {-sz,-sz,h0},  { sz,-sz,h0},  {-sz,sz,h0},   {sz,sz,h0}   }   );
	//terrain_shape   = FieldPatch::makeList( 15, { 0.5,   Vec3d(-tersz,-tersz,h0),  Vec3d( tersz,-tersz,h0),  Vec3d(-tersz,tersz,h0),   Vec3d(tersz,tersz,h0)   }   );

}

/*
int makeSinTerrain( float nx, float ny, float szx, float szy, float height ){
	for (int ix=0; ix<nx; ix++){
		for (int iy=0; iy<nx; iy++){
			val = css[ix][iy];
			glEndList();
	}}
}
*/

void GameWorld::init( ){

	makeEnvironment( 2000.0f );
	printf( " Environment DONE! \n" );
	makeAeroCraft();
    printf( " AeroCraft DONE! \n" );
};




