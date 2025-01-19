
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "globals.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

//int verbosity = 0;

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Vec3Utils.h"

#include "testUtils.h"

#include "SpaceCraft.h"
#include "SpaceCraft2Mesh2.h"
//#include "SoftBody.h"
//#include "Truss.h"
//#include "SpaceCraft2Truss.h" // deprecated

#include "SpaceCraftDraw.h"

//#include "SphereSampling.h"
//#include "DrawSphereMap.h"
//#include "Draw3D_Surf.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"

#include "IO_utils.h"

#include "EditSpaceCraft.h"

#include "TrussDynamics_d.h"
#include "TrussDynamics_f.h"
#include "OCL_Orb.h"
//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

//#include "Tree.h"
//#include "spaceCraftEditorUtils.h"

#include "SpaceCraftGUI.h"
#include "argparse.h"

#include "testUtils.h"

#include "spaceCraftSimulatorOCL.h"
#include "SpaceCraftDynamicsApp.h"

// ======================  Global Variables & Declarations

using namespace SpaceCrafting;

// Mesh::Builder2 mesh;
// OrbSim         sim2;   // CPU double precision
// OCL_Orb        sim_cl; // single precision both CPU and GPU

// // make list of fixed points containing following indices : [0]
// std::vector<int> fixPoints{ 1 };
// double Kfix = 1e+12;

// ======================  Free Functions

// void to_OCL_Orb( OCL_Orb& sim, Mesh::Builder2& mesh, Vec3f p0=Vec3fZero, Vec3f omega=Vec3fZero, bool bCholesky=false, bool bGPU=true ){
//     exportSim( sim, mesh, workshop );

//     sim.reallocFixed();
//     for(int i:fixPoints){ printf("to_OCL_Orb fixing point %i \n", i); sim.kFix[i] = Kfix ; }

//     sim.damping = 1e-5; 
//     sim.set_time_step( 1e-3 );

//     sim.cleanForce();
//     sim.cleanVel();
//     sim.cleanImpuls();
//     sim.addAngularVelocity(  p0, omega );

//     for(int i=0; i<sim.nPoint; i++){
//     //     Vec3f r = sim_cl.points[i].f - p0;
//     //     sim.vel[i].f.set_cross(omega,r);
//            sim.points[i].x -= 100.0;
//     }
//     if(bCholesky){
//         //double omega = 1.0;
//         //double dt    = 0.05;
//         //double dt    = 0.02;
//         //double dt    = 0.01;
//         int n = sim_cl.nPoint;
//         mat2file<float>( "bond_params_f.log",   sim.nBonds,4, (float*) sim.bparams        );
//         mat2file<int>( "neighs_before.log",  n, sim.nNeighMax, sim.neighs, "%5i " );
//         sim.prepare_LinearSystem( sim.dt, true, true, true, 32 );
//         mat2file<int>  ( "neighs_after_f.log", n, sim.nNeighMaxLDLT, sim.neighsLDLT, "%5i " );
//         mat2file<float>( "PDmat_f.log",  n,n, sim.PDmat  );
//         mat2file<float>( "LDLT_L_f.log", n,n, sim.LDLT_L );
//         mat2file<float>( "LDLT_D_f.log", n,1, sim.LDLT_D );
//         //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CholeskySparse;
//         sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Cholesky;
//         //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CG;
//     }
//     if(bGPU){
//         printf("###### OpenCL initialization\n");
//         sim.makeKrenels_Orb( "./common_resources/cl" );
//         sim.initCLBuffsOrb(  );
//         //sim_cl.setup_test_enque();
//         //sim_cl.setup_blur();
//         //sim_cl.test_enque();
//         sim.setKngs();
//         sim.upload( sim_cl.ibuff_points, sim_cl.points );
//         sim.upload( sim_cl.ibuff_vels  , sim_cl.vel    );
//         sim.upload( sim_cl.ibuff_forces, sim_cl.forces );
//         sim.upload( sim_cl.ibuff_impuls, sim_cl.impuls );
//         sim.upload( sim_cl.ibuff_kngs  , sim_cl.kngs   );
//         sim.setup_evalTrussForce1();
//         sim.setup_evalTrussForce2();
//         sim.setup_move();
//         sim.setup_assembleAndMove();
//         sim.setup_evalTrussBondForce();
//     }
// }

// ====================== Class Definitions

class SpaceCraftDynamicsOCLapp : public SpaceCraftDynamicsApp { public:
    //Mesh::Builder2* _mesh=0;
    //OrbSim*         _sim =0;
    OCL_Orb*          _sim_cl=0;
    PickerUI          picker;
    //bool bDouble = true;
    bool bGPU    = true;
    // ==== function declarations

    //void drawSim( OrbSim& sim );
    //void drawSim_f( OrbSim_f& sim );
    //void drawSim( OCL_Orb& sim );
    //void initSim( Mesh::Builder2& mesh );

    //virtual void initSimDefault() override;
	virtual void draw       () override;
	SpaceCraftDynamicsOCLapp( int& id, int WIDTH_, int HEIGHT_ );

    virtual void bindSimulators( SpaceCraftSimulator* simulator_ ){
        simulator=simulator_;
        _mesh  = simulator->getMesh();
        _sim   = simulator->getOrbSim();
        _sim_f = simulator->getOrbSim_f();
        _sim_cl= &((SpaceCraftSimulatorOCL*)simulator)->sim_cl;
    }

};

// void SpaceCraftDynamicsOCLapp::initSimDefault( ){
//     init_workshop ();
//     makeTrussShape( *_mesh, 3, 10, 100.0, 10.0 );
//     initSim( *_mesh );
// }

// void SpaceCraftDynamicsOCLapp::initSim( Mesh::Builder2& mesh ){
//     to_OrbSim  ( *_sim,  mesh );
//     to_OrbSim_f( *_sim_cl, mesh );
//     //to_OCL_Orb( *_sim_cl, mesh );
// }


// void SpaceCraftDynamicsOCLapp::drawSim( & sim_f ){
//     OCL_Orb& sim = *_sim_cl;
//     if( (sim.nPoint==0)||(sim.points==0) ){ printf( "ERROR in runSim() sim.nPoint=%i sim.points=%p => exit() \n", sim.nPoint, sim.points ); exit(0); };
//     //niter=1;
//     int nbig = 20;
//     int nsub = 20;
//     if(bRun){
//         long t0 = getCPUticks();
//         if(bGPU){
//             sim.run_PDcl( nbig, nsub, 0b111, 0b111 );
//         }else{
//             sim.run_Cholesky_omp_simd(nbig);
//         }
//         double T = (getCPUticks()-t0)*1e-6;
//         //printf( "TIME run: %g [Mticks] LinSolve: %g(%g\%) Mdot: %g(%g\%) niter_av=%g err2tot=%g \n", time_run_LinSolve, sim.time_LinSolver,   100*sim.time_LinSolver/time_run_LinSolve, sim.time_cg_dot, 100*sim.time_cg_dot/time_run_LinSolve, sim.cgSolver_niterdone/((double)perFrame), sim.cgSolver.err2tot );
//     }
//     sim.evalBondTension();
//     glColor3f(1.0,0.0,1.0);
//     renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-1 );
//     glColor3f(0.0,1.0,1.0);
//     renderPointForces( sim.nPoint, sim.points, sim.vel, 1e-1 );
//     renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
//     if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
// }

SpaceCraftDynamicsOCLapp::SpaceCraftDynamicsOCLapp( int& id, int WIDTH_, int HEIGHT_ ) : SpaceCraftDynamicsApp( id, WIDTH_, HEIGHT_ ) {  }

void SpaceCraftDynamicsOCLapp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    //glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    //drawBody();
    if(bRun){
        long t0 = getCPUticks(); 
        if( bDouble ){
            _sim->run_Cholesky_omp_simd(perFrame);
            drawSim( *_sim   );
        }else if (bGPU){
            int nsub = 10;
            _sim_cl->run_PDcl( perFrame, nsub, 0b111, 0b111 );
        }else{
            _sim_f->run_Cholesky_omp_simd(perFrame);
        }
        drawSim_f( *_sim_f );
        double T = (getCPUticks()-t0);
        printf( "SpaceCraftDynamicsApp::drawSim(bDoublke=%i) perFrame: %3i nPoint:%6i TIME: %8.3f [Mticks] %8.1f [tick/point] \n", bDouble, perFrame, _sim->nPoint, T*1e-6,  T/(perFrame*_sim->nPoint) );
    }
	//if(!bDrawTrj)glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	//glEnable(GL_DEPTH_TEST);
    Draw3D::drawAxis( 1.0 );
};

// void SpaceCraftDynamicsOCLapp::drawHUD(){
//     glDisable( GL_LIGHTING );
//     glDisable(GL_DEPTH_TEST);
//     //gui.draw();
//     //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
//     //glPopMatrix();
// }

// ===================== MAIN

SpaceCraftSimulatorOCL W;

int main(int argc, char *argv[]){

    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	SpaceCraftDynamicsOCLapp * app = new SpaceCraftDynamicsOCLapp( junk, dm.w-150, dm.h-100 );
    app->bindSimulators( &W ); 

    LambdaDict funcs;
    funcs["-s"]={1,[&](const char** ss){ 
        W.reloadShip( ss[0] );
        W.initSimulators( );
    }};
    funcs["-float" ]={0,[&](const char** ss){ app->bDouble=false; } };
    funcs["-double"]={0,[&](const char** ss){ app->bDouble=true;  } };
    funcs["-gpu"   ]={0,[&](const char** ss){ app->bGPU   =true;  } };
    funcs["-cpu"   ]={0,[&](const char** ss){ app->bGPU   =false;  } };
    funcs["-fix"   ]={1,[&](const char** ss){ readlist( ss[0], W.fixPoints);  } };
    process_args( argc, argv, funcs );
    if( W.sim.nPoint == 0 ){ W.initSimDefault(); }
	app->loop( 1000000 );
	return 0;
}
















