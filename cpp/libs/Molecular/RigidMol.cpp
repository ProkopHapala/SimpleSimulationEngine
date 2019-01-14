

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
    world.gridFF.grid.printCell();
    world.gridFF.loadXYZ     ( fname, params );
    world.translate          ( *(Vec3d*)pos0 );
    world.genPLQ();
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

void debugSaveGridFF(const char* fname, double* testREQ ){
    Vec3d * FFtot = new Vec3d[world.gridFF.grid.getNtot()];
    world.gridFF.evalCombindGridFF( *(Vec3d*)testREQ, FFtot );
    saveXSF( "FFtot_z.xsf", world.gridFF.grid, FFtot, 2, world.gridFF.natoms, world.gridFF.apos, world.gridFF.atypes );
    delete [] FFtot;
}

// ========= Molecule initialization

void initParams( char* fname_atomTypes, char* fname_bondTypes ){
    builder.params = &params;
    if(fname_atomTypes) params.loadAtomTypes( fname_atomTypes );
    if(fname_bondTypes) params.loadBondTypes( fname_bondTypes );
    printf("initParams done! \n");
}

int loadMolType   ( const char* fname ){ return builder.loadMolType(fname ); };
int insertMolecule( int itype, double* pos, double* rot, bool rigid ){ return builder.insertMolecule( itype, *(Vec3d*)pos, *(Mat3d*)rot, rigid ); };

void bakeMMFF(){
    world.printAtomInfo();
    builder.toMMFF( &world );                                 DEBUG
    world.printAtomInfo(); //exit(0);
    //world.allocFragment( nFrag );
    //opt.bindArrays( 8*world.nFrag, (double*)world.poses, new double[8*world.nFrag], (double*)world.poseFs ); 
}

void prepareOpt(){
    //opt.bindArrays( 8*world.nFrag, world.poses, world.poseVs, world.poseFs );
    world.allocateDyn(); 
    world.initDyn();     
    opt.bindArrays( world.nDyn, world.dynPos, world.dynVel, world.dynForce ); DEBUG
    opt.setInvMass( 1.0 );
    opt.cleanVel  ( );
    //exit(0);
    //printf("POSE_pos   : \n"); printPoses( world.nFrag, world.poses  );
    //printf("POSE_Force : \n"); printPoses( world.nFrag, world.poseFs );
}

double relaxNsteps( int nsteps, double F2conf ){
    double F2=1e+300;
    for(int itr=0; itr<nsteps; itr++){
        world.cleanAtomForce();
        world.frags2atoms();
        world.eval_FFgrid();
        world.eval_MorseQ_On2_fragAware();

        world.cleanPoseTemps();
        world.aforce2frags();   

        //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); } // gradient descent
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        F2 = opt.move_FIRE();
        if(F2<F2conf) break;

        //world.toDym(true);
        //world.checkPoseUnitary();
        //world.fromDym();
    }
    return F2;
}

void save2xyz( char * fname ){
    save2xyz( fname, &world, &params ); 
}

} // extern "C"{
