

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
//#include "Draw.h"
//#include "Draw3D.h"
//#include "SDL_utils.h"
//#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "raytrace.h"
#include "Molecule.h"
#include "MMFF.h"
#include "MMFFBuilder.h"

#include "IO_utils.h"
#include "geom3D.h"

#include "DynamicOpt.h"

#include "AtomicConfiguration.h"
#include "DistanceHierarchy.h"

//#include "AppSDL2OGL_3D.h"
//#include "MolecularDraw.h"


#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Multipoles.h"
#include "DynamicOpt.h"

#include "radial_splines.h"
#include "AtomTypes.h"
#include "MoleculeType.h"
#include "MolecularWorld.h"

#include "PointCloudComparator.h"
#include "SphereTreeND.h"

#include "FitFF.h"

//#include "testUtils.h"

// ============ Global Variables

//std::vector<Molecule*> molTypes;
MMFFparams  params;
MMFF        world;
MMFFBuilder builder;

DynamicOpt  opt;

std::vector<FILE*> files;

//FastAtomicMetric    atomdist;
//AtomicConfiguration conf1;
//DistanceHierarchy<AtomicConfiguration> database;
//AtomicManipulator manipulator;

extern "C"{

// ========= Grid initialization

void initRigidSubstrate( char* fname, int* ns, double* pos0, double* cell ){
    //printf( "params.atypNames:\n" );
    //for(auto kv : params.atypNames) { printf(" %s %i \n", kv.first.c_str(), kv.second ); }
    world.gridFF.grid.n    = *(Vec3i*)ns;
    world.gridFF.grid.pos0 = *(Vec3d*)pos0;
    //world.gridFF.loadCell ( "inputs/cel.lvs" );
    world.gridFF.grid.setCell( *(Mat3d*)cell );
    //world.gridFF.grid.printCell();
    world.gridFF.loadXYZ     ( fname, params );
    world.translate          ( *(Vec3d*)pos0 );
    world.gridFF.allocateFFs();
    //world.gridFF.evalGridFFs( {0,0,0} );
}

void recalcGridFF( int* ns){
    world.gridFF.evalGridFFs( *(Vec3i*)ns );
}

void saveGridFF(){
    if(world.gridFF.FFelec )  saveBin( "data/FFelec.bin",   world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFelec );
    if(world.gridFF.FFPauli)  saveBin( "data/FFPauli.bin",  world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFPauli );
    if(world.gridFF.FFLondon) saveBin( "data/FFLondon.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFLondon );
}

void loadGridFF(){
    if(world.gridFF.FFelec )  loadBin( "data/FFelec.bin",   world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFelec );
    if(world.gridFF.FFPauli)  loadBin( "data/FFPauli.bin",  world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFPauli );
    if(world.gridFF.FFLondon) loadBin( "data/FFLondon.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFLondon );
}

void debugSaveGridFF( char* fname, double* testREQ ){
    Vec3d * FFtot = new Vec3d[world.gridFF.grid.getNtot()];
    world.gridFF.evalCombindGridFF( *(Vec3d*)testREQ, FFtot );
    saveXSF( fname, world.gridFF.grid, FFtot, 2, world.gridFF.natoms, world.gridFF.apos, world.gridFF.atypes );
    delete [] FFtot;
}

// ========= Molecule initialization

void initParams( char* fname_atomTypes, char* fname_bondTypes ){
    builder.params = &params;
    if(fname_atomTypes) params.loadAtomTypes( fname_atomTypes );
    if(fname_bondTypes) params.loadBondTypes( fname_bondTypes );
    printf( "params.atypNames.size() %i \n", params.atypNames.size() );
    //printf("initParams done! \n");
}

int loadMolType   ( char* fname ){ return builder.loadMolType(fname ); };
int insertMolecule( int itype, double* pos, double* rot, bool rigid ){ return builder.insertMolecule( itype, *(Vec3d*)pos, *(Mat3d*)rot, rigid ); };

void clear(){
    builder.clear();
}

void bakeMMFF(){
    builder.toMMFF( &world );
    world.genPLQ();
    world.printAtomInfo(); //exit(0);
    //world.allocFragment( nFrag );
    //opt.bindArrays( 8*world.nFrag, (double*)world.poses, new double[8*world.nFrag], (double*)world.poseFs ); 
}

void prepareOpt(){
    //opt.bindArrays( 8*world.nFrag, world.poses, world.poseVs, world.poseFs );
    printf("DEBUG a.0\n");
    world.allocateDyn();   printf("DEBUG a.1\n");
    world.initDyn();       printf("DEBUG a.2\n");
    opt.bindArrays( world.nDyn, world.dynPos, world.dynVel, world.dynForce );  printf("DEBUG a.3\n");
    opt.setInvMass( 1.0 ); printf("DEBUG a.4\n");
    opt.cleanVel  ( );     printf("DEBUG a.5\n");
    //exit(0);
    //printf("POSE_pos   : \n"); printPoses( world.nFrag, world.poses  );
    //printf("POSE_Force : \n"); printPoses( world.nFrag, world.poseFs );
    DEBUG
}

void setOptFIRE(
	double dt_max,
	double dt_min,
	double damp_max,
	int    minLastNeg,
	double finc,
	double fdec,
	double falpha,
	double kickStart
    ){
	opt.minLastNeg   = minLastNeg;
	opt.finc         = finc;
	opt.fdec         = fdec;
	opt.falpha       = falpha;
	opt.kickStart    = kickStart;
	opt.dt_max       = dt_max;
	opt.dt_min       = dt_min;
	opt.damp_max     = damp_max;
}

double relaxNsteps( int nsteps, double F2conf ){
    double F2=1e+300;
    //DEBUG
    for(int itr=0; itr<nsteps; itr++){
        //printf( "===== relaxNsteps itr %i \n", itr );
        world.cleanAtomForce();
        world.frags2atoms();
        if( world.gridFF.FFPauli ) world.eval_FFgrid();
        world.eval_MorseQ_On2_fragAware();

        world.cleanPoseTemps();
        world.aforce2frags();

        world.toDym(true);
        //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); } // gradient descent
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        F2 = opt.move_FIRE();
        //printf( "F2 %g dt %g \n", F2, opt.dt );
        if(F2<F2conf) break;
        world.checkPoseUnitary();
        world.fromDym();
        //DEBUG

        //printf( ">> itr %i F2 %g dt %g qrot (%g,%g,%g,%g) int %li \n", itr, F2, opt.dt, world.poses[4], world.poses[5], world.poses[6], world.poses[7], world.gridFF.FFPauli );
        printf( ">> itr %i F2 %g dt %g poses (%g,%g,%g,%g, %g,%g,%g,%g) \n", itr, F2, world.poses[0], world.poses[1], world.poses[2], world.poses[3], world.poses[4], world.poses[5], world.poses[6], world.poses[7] );

    }
    return F2;
}

void save2xyz( char * fname ){
    save2xyz( fname, &world, &params ); 
}

void write2xyz( int i ){
     if( (i>=0)&&(i<files.size()) ){ write2xyz( files[i], &world, &params );  }
}

void closef(int i){ if( (i>=0)&&(i<files.size()) ) if( files[i] ) fclose( files[i] ); }

int openf(char* fname, int i, char* mode ){
    if((i<0)||(i>=files.size())){
        i = files.size();
        files.push_back(NULL);
    }else{
        closef(i);
    }
    files[i] = fopen(fname, "w");
    return files.size()-1;
};

double* getPoses  (int* n){ *n=world.nFrag;  return world.poses;         }
double* getAtomPos(int* n){ *n=world.natoms; return (double*)world.apos; }

} // extern "C"{
