
#include <SDL2/SDL_opengl.h>

#include "Draw3D.h"

#include "testUtils.h"

#include "AeroCraftWorld.h" // THE HEADER
#include "testUtils.h"

void AeroCraftWorld::reallocateTrj(int n){
  ntrjMax=n;
  if(trjPos  ) delete trjPos;   trjPos   = new Vec3d[n];
  if(trjVel  ) delete trjVel;   trjVel   = new Vec3d[n];
  if(trjForce) delete trjForce; trjForce = new Vec3d[n];
  if(trjFw   ) delete trjFw;    trjFw    = new Vec3d[n];
  if(trjUp   ) delete trjUp;    trjUp    = new Vec3d[n];
  if(trjT    ) delete trjT;     trjT     = new double[n];
}



void AeroCraftWorld::resetSteer( ){
    myCraft->panels[0].lrot = myCraft_bak->panels[0].lrot;
    myCraft->panels[1].lrot = myCraft_bak->panels[1].lrot;
    myCraft->panels[2].lrot = myCraft_bak->panels[2].lrot;
    myCraft->panels[3].lrot = myCraft_bak->panels[3].lrot;
}

void AeroCraftWorld::steerToDir( const Vec3d& dir ){

    Mat3d rotMatT;
    rotMatT.setT(myCraft->rotMat);
    Draw3D::drawMatInPos( rotMatT, myCraft->pos );
    double dyaw   = rotMatT.a.dot( dir );
    double dpitch = rotMatT.b.dot( dir );
    const double acut = 0.1;
    //double droll = (a>acut)?(a-acut):((a<-acut)?(a+acut):0.0d);
    double droll = dyaw;

    /*
	glColor4f(1.0f,1.0f,1.0f,0.9f);
	char str[256];
	sprintf(str, "a %3.3f b 3.3f %3.3f\0",a,b);
	Draw3D::drawText(str, myCraft->pos, fontTex_DEBUG, 0.2, 0, 0 );
	glEnable(GL_DEPTH_TEST);
	*/
    resetSteer();
    myCraft->steerTo( 0.1*droll, 0.5*dpitch+0.5*fabs(dyaw), 0.5*dyaw);
    /*
    Draw3D::drawMatInPos( myCraft->rotMat, myCraft->pos );
    double a = myCraft->rotMat.a.dot( dir );
    double b = myCraft->rotMat.b.dot( dir );
    myCraft->panels[0].lrot.rotate(  _clamp(  0.2*a             , 0.0, 0.5),  {1.0,0.0,0.0} );
    myCraft->panels[1].lrot.rotate(  _clamp( -0.2*a             , 0.0, 0.5),  {1.0,0.0,0.0} );
    myCraft->panels[2].lrot.rotate(  _clamp(  1.5*fabs(a)-1.5*b ,-0.5, 0.5), {1.0,0.0,0.0} );
    myCraft->panels[3].lrot.rotate(  _clamp(  1.5*a             ,-0.5, 0.5),   {0.0,1.0,0.0} );
    */
};

/*
void AeroCraftWorld::steerToMat( Mat3d& mat, bool on ){
    if(on) resetSteer();
    Mat3d rotMatT;
    rotMatT.setT(myCraft->rotMat);
    Draw3D::drawMatInPos( rotMatT, myCraft->pos );
    double a = rotMatT.a.dot( dir );
    double b = rotMatT.b.dot( dir );
    const double acut = 0.1;
    double a_ = (a>acut)?(a-acut):((a<-acut)?(a+acut):0.0d);

	glColor4f(1.0f,1.0f,1.0f,0.9f);
	char str[256];
	sprintf(str, "a %3.3f b 3.3f %3.3f\0",a,b);
	Draw3D::drawText(str, myCraft->pos, fontTex_DEBUG, 0.2, 0, 0 );
	glEnable(GL_DEPTH_TEST);

	if(on){

        myCraft->panels[0].lrot.rotate(  _clamp( -0.05*a_            , 0.0, 0.5),  {1.0,0.0,0.0} );
        myCraft->panels[1].lrot.rotate(  _clamp( +0.05*a_            , 0.0, 0.5),  {1.0,0.0,0.0} );
        myCraft->panels[2].lrot.rotate(  _clamp(  0.5*b + 0.5*fabs(a),-0.5, 0.5), {1.0,0.0,0.0} );
        myCraft->panels[3].lrot.rotate(  _clamp( -0.5*a ,-0.5, 0.5),   {0.0,1.0,0.0} );
    }
};
*/

void AeroCraftWorld::update( ){

	//long long nt1 = getCPUticks();


	myCraft->checkStateNormal();
	//for(int iSubStep=0; iSubStep<PHYS_STEPS_PER_FRAME; iSubStep++ ){
	for( int i=0; i<perFrame; i++ ){
		myCraft->clean_temp();
		myCraft->force.set    ( { 0, gravityG*myCraft->mass, 0 } );
		//myCraft->force.add_mul(   myCraft->rotMat.c, 500.0 ); // motor
		myCraft->applyAeroForces( {0,0,0} );
		myCraft->move(dt);
	}

	//long long nt2  = getCPUticks();
	//double perStep = double(nt2-nt1)/PHYS_STEPS_PER_FRAME;
	//tickSum+= perStep; stepSum++;
	//printf( " PERFORMANCE: %f ticks/step ( in %i steps ) average: %f \n", perStep, PHYS_STEPS_PER_FRAME, tickSum/stepSum );

};

/*
void AeroCraftWorld::makeAeroCraft(){
	const int len = 5;
	//                      motor        wingLeft   wingRight  elevator   rudder
	double masses[len]  = { 4,           1,         1,         0.5f,      0.5f      };
	Vec3d  poss[len]    = { {0,-0.1,1}, {-2,0,0}, {2,0,0},   {0,0,-3},   {0,0.2,-3} };
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

	myCraft->wingLeft .lrot.set( { 1,0,0, 0,1,0,  0,0,1  } ); myCraft->wingLeft .lrot.rotate( -0.1, { 0,0,1 } );
	myCraft->wingRight.lrot.set( { 1,0,0, 0,1,0,  0,0,1  } ); myCraft->wingRight.lrot.rotate( +0.1, { 0,0,1 } );
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
*/


int AeroCraftWorld::makeBuildingsGrid( int nx, int ny, float sx, float sy, float cx, float cy,  float min_height, float max_height ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
	for (int ix=-nx; ix<nx; ix++){
		float x = ix*sx;
		for (int iy=-ny; iy<ny; iy++){
			float height = randf() * (max_height-min_height) + min_height;
			float y = iy*sy;
			Draw3D::drawBox( x, x + sx*cx, 0, height, y, y + sy*cy, 0.75f, 0.75f, 0.75f );
		}
	}
	glEndList();
	return( ilist );
}

int AeroCraftWorld::makeBuildingsClusters( int nclustest, int nmin, int nmax, float minx, float maxx, float miny, float maxy, float min_dist, float max_dist, float min_size, float max_size, float min_height, float max_height ){
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
			Draw3D::drawBox( x-dx, x+dx, 0, height, y-dy, y+dy, 0.75f, 0.75f, 0.75f );
			nboxes++;
		};
	};
	printf(" %i buildings \n", nboxes );
	glEndList();
	return( ilist );
}

void AeroCraftWorld::makeEnvironment( float sz ){

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
	terrain_shape     = fieldPatch.makeList( 5, { 0.5d,   {-sz,-sz,h0},  { sz,-sz,h0},  {-sz,sz,h0},   {sz,sz,h0}   }   );
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

void AeroCraftWorld::evalAircraftTrajectory( int n, int nsub, int msub, double dt ){
    ntrj=n;
    double t=0;
    long ticks1 = getCPUticks();
    for(int i=0; i<n; i++){
        trjPos  [i] = myCraft->pos;
        trjVel  [i] = myCraft->vel;
        trjForce[i] = myCraft->force;
        trjFw   [i] = myCraft->rotMat.c;
        trjUp   [i] = myCraft->rotMat.b;
        trjT    [i] = t;
        for(int j=0; j<nsub; j++){
            //resetSteer();
            autoPilot1.control(dt*msub);
            for(int itr=0; itr<msub; itr++){
                myCraft->clean_temp();
                myCraft->force.set      ( { 0, gravityG*myCraft->mass, 0 } );
                myCraft->applyAeroForces( {0,0,0} );
                myCraft->move(dt);
                t+=dt;
            }
		}
    }
    double  dticks        = getCPUticks() - ticks1;
    double  ticksPerIter  = dticks/(n*nsub*msub);
    double  tickPerSec    = dticks/(t-trjT[0]);
    printf( "Ticks %g /iter %g /sec %g\n", dticks, ticksPerIter, tickPerSec );
}

void AeroCraftWorld::doStaticTesting( ){
    //--- propeller characterisic
    double vmin=0.01;
    double vmax=300.0;
    double dv  =5.0;
    for(double v=vmin; v<vmax; v+=dv){
        double thrust = myCraft->propelers[0].getThrust(v);
        printf(" v=%f [m/s] thrust=%f [N] \n",  v, thrust );
    }

    int ntrj_=500;
    reallocateTrj(ntrj_);
    evalAircraftTrajectory(ntrj_, 100, 10, 0.001);  // correspond to interactive simulation
    //evalAircraftTrajectory(ntrj_, 25, 4, 0.01);
    //exit(0);
}

void AeroCraftWorld::init( ){

	makeEnvironment( 20000.0f );
	printf( " Environment DONE! \n" );
	//makeAeroCraft();

	myCraft_bak = new AeroCraft();   myCraft_bak->fromFile("data/AeroCraft1.ini");
    myCraft     = new AeroCraft();   myCraft    ->fromFile("data/AeroCraft1.ini");

    myCraft->pos.y=200.0;

    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );

    autoPilot1.craft=myCraft;

    staticTest=false;

    // static test ---------
    if( staticTest ) doStaticTesting( );

    printf( " AeroCraft DONE! \n" );
};




