#ifndef SpaceCraftSimulatorOCL_h
#define SpaceCraftSimulatorOCL_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "SpaceCraft.h"
#include "SpaceCraft2Mesh2.h"
#include "SpaceCraftDraw.h"
#include "OrbSim_d.h"
#include "OrbSim_f.h"
#include "OCL_Orb.h"
#include "spaceCraftSimulator.h"

namespace SpaceCrafting {

void distort_points( int n, Quat4f* ps, Quat4f* ps_out, float rnd=0.1, Vec3f sc=Vec3fOne, int seed=15454 ){
    srand(seed);
    for(int i=0; i<n; i++){
        Vec3f r; r.fromRandomCube(rnd);
        ps_out[i].f = ps[i].f*sc + r;
    }
}

void test_SolverConvergence( OCL_Orb& sim_cl, Quat4f* ps_bak , int nSolverIters=10, int nbmix=3, float bmix0=0.5, float dbmix=0.1 ){
    printf( "test_SolverConvergence() nPoint=%i nSolverIters=%g nbmix=%i bmix0=%g dbmix=%g \n", sim_cl.nPoint, nSolverIters, nbmix, bmix0, dbmix );
    if( ps_bak){ for(int i=0; i<sim_cl.nPoint; i++){ sim_cl.points[i]=ps_bak[i]; } }
    sim_cl.run_SolverConvergence( nSolverIters, 0, true );
    if( ps_bak){ for(int i=0; i<sim_cl.nPoint; i++){ sim_cl.points[i]=ps_bak[i]; } }
    sim_cl.run_SolverConvergence( nSolverIters, 2, true );
    for(int i=0; i<nbmix; i++){
        sim_cl.bmix.y = bmix0 + i*dbmix; 
        if( ps_bak){ for(int i=0; i<sim_cl.nPoint; i++){ sim_cl.points[i]=ps_bak[i]; } }
        sim_cl.run_SolverConvergence( nSolverIters, 1, true );
    }
}

void initGPU( OCL_Orb& sim ){
    printf("#### ====  SpaceCraftSimulatorOCL::to_OCL_Orb() \n");
    if( sim.nPoint==0 ){ printf( "ERROR in SpaceCraftSimulatorOCL::to_OCL_Orb() sim.nPoint=%i => run  SpaceCraftSimulator::to_OrbSim_f() first \n... exit() \n", sim.nPoint ); exit(0); };
    if( !sim.bOpenCL_initialized ){
        sim.print_devices();
        sim.initOCL();
    }
    sim.cleanImpuls();
    printf("###### OpenCL initialization\n");
    sim.makeKrenels_Orb( "./common_resources/cl" );
    sim.initCLBuffsOrb(  );
    //sim_cl.setup_test_enque();
    //sim_cl.setup_blur();
    //sim_cl.test_enque();
    sim.setKngs();
    sim.upload( sim.ibuff_points, sim.points );
    sim.upload( sim.ibuff_vels  , sim.vel    );
    sim.upload( sim.ibuff_forces, sim.forces );
    sim.upload( sim.ibuff_impuls, sim.impuls );
    sim.upload( sim.ibuff_kngs  , sim.kngs   );
    sim.setup_evalTrussForce1();
    sim.setup_evalTrussForce2();
    sim.setup_move();
    sim.setup_assembleAndMove();
    sim.setup_evalTrussBondForce();
}

// void init_workshop(){
//     //StickMaterial *o = new StickMaterial();
//     workshop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
//     workshop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );
// }

// =======================================================================
class SpaceCraftSimulatorOCL: public SpaceCraftSimulator { public:

    // Mesh::Builder2 mesh;
    // OrbSim         sim;
    OCL_Orb        sim_cl;
    // std::vector<int> fixPoints;
    // bool   bDouble = false;
    // bool   bRun    = false;
    // double time=0;

    // ---- Functions
    // void         reloadShip ( const char* fname );
    // virtual void to_OrbSim  ( double dt=0.1, Vec3d p0=Vec3dZero, Vec3d omega=Vec3dZero );
    // virtual void to_OrbSim_f( double dt=0.1, Vec3f p0=Vec3fZero, Vec3f omega=Vec3fZero );
    void to_OCL_Orb( Vec3f p0=Vec3fZero, Vec3f omega=Vec3fZero, bool bCholesky=false, bool bGPU=true );
    virtual void initSimulators( double dt=0.1, Vec3d p0=Vec3dZero, Vec3d omega=Vec3dZero ) override;
    // virtual void initSimDefault();

    // virtual OrbSim*         getOrbSim  (){ return &sim;   };
    virtual OrbSim_f*       getOrbSim_f(){ return &sim_cl; };
    // virtual Mesh::Builder2* getMesh(){ return &mesh; };

};
    
void SpaceCraftSimulatorOCL::initSimulators( double dt, Vec3d p0, Vec3d omega ){
    printf("SpaceCraftSimulatorOCL::initSimulators() \n");
    to_OrbSim  ( sim,    mesh, dt,        p0,        omega, fixPoints.size(), fixPoints.data() );
    to_OrbSim_f( sim_cl, mesh, dt, (Vec3f)p0, (Vec3f)omega, fixPoints.size(), fixPoints.data() );
    initGPU    ( sim_cl );
}

}// namespace SpaceCrafting

#endif // SpaceCraftSimulatorOCL_h
