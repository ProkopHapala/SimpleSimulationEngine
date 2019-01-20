
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

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

// ==== Global Variables

MMFFparams  params;
MMFF        world;
MMFFBuilder builder;
DynamicOpt  opt;

//int     fontTex;
//int     ogl_sph;

char str[256];

void   initRigidSubstrate();
void   bakeMMFF();
void   prepareOpt();
double relaxNsteps( int nsteps, double F2conf );

// ==== Function Implementation

void initRigidSubstrate(){
    // ---- Rigid Substrate
    printf( "params.atypNames:\n" );
    for(auto kv : params.atypNames) { printf(" %s %i \n", kv.first.c_str(), kv.second ); }
    world.gridFF.grid.n    = (Vec3i){60,60,100};
    world.gridFF.grid.pos0 = (Vec3d){0.0d,0.0d,0.0d};
    world.gridFF.loadCell ( "inputs/cel.lvs" );
    world.gridFF.grid.printCell();
    world.gridFF.loadXYZ  ( "inputs/NaCl_wo4.xyz", params );
    world.translate( {0.0,0.0,4.5} );

    Vec3d testREQ,testPLQ;
    testREQ = (Vec3d){ 1.487, sqrt(0.0006808), 0.0}; // H
    testPLQ = REQ2PLQ( testREQ, -1.6 );//
    world.genPLQ();
    world.gridFF.allocateFFs();
    bool recalcFF = false;
    if( recalcFF ){
        world.gridFF.evalGridFFs( {1,1,1} );
        if(world.gridFF.FFelec )  saveBin( "data/FFelec.bin",   world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFelec );
        if(world.gridFF.FFPauli)  saveBin( "data/FFPauli.bin",  world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFPauli );
        if(world.gridFF.FFLondon) saveBin( "data/FFLondon.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFLondon );
    }else{
        if(world.gridFF.FFelec )  loadBin( "data/FFelec.bin",   world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFelec );
        if(world.gridFF.FFPauli)  loadBin( "data/FFPauli.bin",  world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFPauli );
        if(world.gridFF.FFLondon) loadBin( "data/FFLondon.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFLondon );
    }
    int iatom = 11;
    printf( "testREQ   (%g,%g,%g) -> PLQ (%g,%g,%g) \n",        testREQ.x, testREQ.y, testREQ.z, testPLQ.x, testPLQ.y, testPLQ.z   );
    printf( "aREQs[%i] (%g,%g,%g) -> PLQ (%g,%g,%g) \n", iatom, world.aREQ[iatom].x, world.aREQ[iatom].y, world.aREQ[iatom].z, world.aPLQ[iatom].x, world.aPLQ[iatom].y, world.aPLQ[iatom].z );
    Vec3d * FFtot = new Vec3d[world.gridFF.grid.getNtot()];
    world.gridFF.evalCombindGridFF( testREQ, FFtot );
    saveXSF( "FFtot_z.xsf", world.gridFF.grid, FFtot, 2, world.gridFF.natoms, world.gridFF.apos, world.gridFF.atypes );
    
    //isoOgl = glGenLists(1);
    //glNewList(isoOgl, GL_COMPILE);
    //    renderSubstrate_( world.gridFF.grid, FFtot, world.gridFF.FFelec, 0.01, true );
    //    Draw3D::drawAxis(1.0);
    //glEndList();
    //cam.pos.z = +5.0;
}

/*
void initParams( char* fname_atomTypes, char* fname_bondTypes ){
    builder.params = &params;
    if(fname_atomTypes) params.loadAtomTypes( fname_atomTypes );
    if(fname_bondTypes) params.loadBondTypes( fname_bondTypes );
    printf( "params.atypNames.size() %i \n", params.atypNames.size() );
    //printf("initParams done! \n");
}
*/

//int loadMolType   ( char* fname ){ return builder.loadMolType(fname ); };
//int insertMolecule( int itype, double* pos, double* rot, bool rigid ){ return builder.insertMolecule( itype, *(Vec3d*)pos, *(Mat3d*)rot, rigid ); };

void bakeMMFF(){
    builder.toMMFF( &world );
    world.genPLQ();
    world.printAtomInfo(); //exit(0);
    //world.allocFragment( nFrag );
    //opt.bindArrays( 8*world.nFrag, (double*)world.poses, new double[8*world.nFrag], (double*)world.poseFs ); 
}

void prepareOpt(){
    //opt.bindArrays( 8*world.nFrag, world.poses, world.poseVs, world.poseFs );
    //printf("DEBUG a.0\n");
    world.allocateDyn();   //printf("DEBUG a.1\n");
    world.initDyn();       //printf("DEBUG a.2\n");
    opt.bindArrays( world.nDyn, world.dynPos, world.dynVel, world.dynForce, world.dynInvMass );  //printf("DEBUG a.3\n");
    opt.setInvMass( 1.0 ); //printf("DEBUG a.4\n");
    opt.cleanVel  ( );     //printf("DEBUG a.5\n");
    //exit(0);
    //printf("POSE_pos   : \n"); printPoses( world.nFrag, world.poses  );
    //printf("POSE_Force : \n"); printPoses( world.nFrag, world.poseFs );
    //DEBUG
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
        F2 = opt.move_FIRE(); 
        //printf( "F2 %g dt %g \n", F2, opt.dt );
        if(F2<F2conf) break;
        world.checkPoseUnitary(); 
        world.fromDym(); 
        //DEBUG
        printf( ">> itr %i F2 %g dt %g qrot (%g,%g,%g,%g) int %li \n", itr, F2, opt.dt, world.poses[4], world.poses[5], world.poses[6], world.poses[7], world.gridFF.FFPauli );
        //printf( ">> itr %i F2 %g dt %g poses (%g,%g,%g,%g, %g,%g,%g,%g) \n", itr, F2, world.poses[0], world.poses[1], world.poses[2], world.poses[3], world.poses[4], world.poses[5], world.poses[6], world.poses[7] );
    }
    return F2;
}

int main(){

    // ======= common Potentials etc.
    builder.params = &params;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );

    Mat3d rot0; rot0.setOne();
    int itype=-1;

    printf( "// =========== System 1 \n" );
    builder.clear();
    itype = builder.loadMolType( "inputs/water_T5_ax.xyz" );
    builder.insertMolecule( itype, (Vec3d){5.78, 6.7, 12.24}, rot0, true );
    bakeMMFF();
    prepareOpt();
    print("DEBUG prepareOpt() -> relaxNsteps \n");
    relaxNsteps( 3, 0.0 );

    // =========== System 2
    printf( "// =========== System 2 \n" );
    builder.clear();
    itype = builder.loadMolType( "inputs/Campher.xyz" );
    builder.insertMolecule( itype, (Vec3d){5.78, 6.7, 12.24}, rot0, true );
    bakeMMFF();
    prepareOpt();
    relaxNsteps( 3, 0.0 );


    for(Molecule* m : builder.molTypes ){ m->dealloc(); delete m; };


    printf( "ALL DONE !\n" );
    exit(0);

}

