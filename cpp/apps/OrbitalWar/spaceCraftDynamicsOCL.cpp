
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

int verbosity = 0;

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

#include "OrbSim_d.h"
#include "OrbSim_f.h"
#include "OCL_Orb.h"
//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

//#include "Tree.h"
//#include "spaceCraftEditorUtils.h"

#include "SpaceCraftGUI.h"
#include "argparse.h"

#include "testUtils.h"

#include "spaceCraftDynamicsShared.h"

// ======================  Global Variables & Declarations

using namespace SpaceCrafting;

Mesh::Builder2 mesh;
OrbSim         sim2;   // CPU double precision
OCL_Orb        sim_cl; // single precision both CPU and GPU

// make list of fixed points containing following indices : [0]
std::vector<int> fixPoints{ 1 };
double Kfix = 1e+12;

// ======================  Free Functions

void to_OCL_Orb( OCL_Orb& sim, Mesh::Builder2& mesh, Vec3f p0=Vec3fZero, Vec3f omega=Vec3fZero, bool bCholesky=false, bool bGPU=true ){
    exportSim( sim, mesh, workshop );

    sim.reallocFixed();
    for(int i:fixPoints){ printf("to_OCL_Orb fixing point %i \n", i); sim.kFix[i] = Kfix ; }

    sim.damping = 1e-5; 
    sim.set_time_step( 1e-3 );

    sim.cleanForce();
    sim.cleanVel();
    sim.cleanImpuls();
    sim.addAngularVelocity(  p0, omega );

    for(int i=0; i<sim.nPoint; i++){
    //     Vec3f r = sim_cl.points[i].f - p0;
    //     sim.vel[i].f.set_cross(omega,r);
           sim.points[i].x -= 100.0;
    }
    if(bCholesky){
        //double omega = 1.0;
        //double dt    = 0.05;
        //double dt    = 0.02;
        //double dt    = 0.01;
        int n = sim_cl.nPoint;
        mat2file<float>( "bond_params_f.log",   sim.nBonds,4, (float*) sim.bparams        );
        mat2file<int>( "neighs_before.log",  n, sim.nNeighMax, sim.neighs, "%5i " );
        sim.prepare_LinearSystem( sim.dt, true, true, true, 32 );
        mat2file<int>  ( "neighs_after_f.log", n, sim.nNeighMaxLDLT, sim.neighsLDLT, "%5i " );
        mat2file<float>( "PDmat_f.log",  n,n, sim.PDmat  );
        mat2file<float>( "LDLT_L_f.log", n,n, sim.LDLT_L );
        mat2file<float>( "LDLT_D_f.log", n,1, sim.LDLT_D );
        //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CholeskySparse;
        sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Cholesky;
        //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CG;
    }
    if(bGPU){
        printf("###### OpenCL initialization\n");
        sim.makeKrenels_Orb( "./common_resources/cl" );
        sim.initCLBuffsOrb(  );
        //sim_cl.setup_test_enque();
        //sim_cl.setup_blur();
        //sim_cl.test_enque();
        sim.setKngs();
        sim.upload( sim_cl.ibuff_points, sim_cl.points );
        sim.upload( sim_cl.ibuff_vels  , sim_cl.vel    );
        sim.upload( sim_cl.ibuff_forces, sim_cl.forces );
        sim.upload( sim_cl.ibuff_impuls, sim_cl.impuls );
        sim.upload( sim_cl.ibuff_kngs  , sim_cl.kngs   );
        sim.setup_evalTrussForce1();
        sim.setup_evalTrussForce2();
        sim.setup_move();
        sim.setup_assembleAndMove();
        sim.setup_evalTrussBondForce();
    }
}

void distort_points( int n, Quat4f* ps, Quat4f* ps_out, float rnd=0.1, Vec3f sc=Vec3fOne, int seed=15454 ){
    srand(seed);
    for(int i=0; i<n; i++){
        Vec3f r; r.fromRandomCube(rnd);
        ps_out[i].f = ps[i].f*sc + r;
    }
}

void test_SolverConvergence( Mesh::Builder2& mesh, Quat4f* ps_bak , int nSolverIters=10, int nbmix=3, float bmix0=0.5, float dbmix=0.1 ){
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

// ====================== Class Definitions

class SpaceCraftDynamicsOCLapp : public SpaceCraftDynamicsApp { public:

    //Mesh::Builder2* _mesh=0;
    //OrbSim*         _sim =0;
    OCL_Orb*          _sim_cl=0;
    PickerUI          picker;

    bool bDouble = true;
    bool bGPU    = true;

    // ==== function declarations

    //void drawSim( OrbSim& sim );
    void drawSim( OCL_Orb& sim );
    void initSim( Mesh::Builder2& mesh );

    virtual void initSimDefault() override;
	virtual void draw          () override;
	//virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	//virtual void eventHandling   ( const SDL_Event& event  ) override;
	//virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );
	SpaceCraftDynamicsOCLapp( int& id, int WIDTH_, int HEIGHT_ );

};

void SpaceCraftDynamicsOCLapp::initSimDefault( ){
    init_workshop ();
    makeTrussShape( *_mesh, 3, 10, 100.0, 10.0 );
    initSim( *_mesh );
}

void SpaceCraftDynamicsOCLapp::initSim( Mesh::Builder2& mesh ){
    to_OrbSim ( *_sim,    mesh );
    to_OCL_Orb( *_sim_cl, mesh );
}


void SpaceCraftDynamicsOCLapp::drawSim( OCL_Orb& sim ){
    if( (sim.nPoint==0)||(sim.points==0) ){ printf( "ERROR in runSim() sim.nPoint=%i sim.points=%p => exit() \n", sim.nPoint, sim.points ); exit(0); };
    //niter=1;
    int nbig = 20;
    int nsub = 20;
    //niter = nbig*nsub;
    sim.set_time_step(0.0001);
    if(bRun){
        long t0 = getCPUticks();
        if(bGPU){
            sim.run_PDcl( nbig, nsub, 0b111, 0b111 );
        }else{
            sim.run_Cholesky_omp_simd(nbig);
        }
        double T = (getCPUticks()-t0)*1e-6;
        //printf( "TIME run: %g [Mticks] LinSolve: %g(%g\%) Mdot: %g(%g\%) niter_av=%g err2tot=%g \n", time_run_LinSolve, sim.time_LinSolver,   100*sim.time_LinSolver/time_run_LinSolve, sim.time_cg_dot, 100*sim.time_cg_dot/time_run_LinSolve, sim.cgSolver_niterdone/((double)perFrame), sim.cgSolver.err2tot );
    }
    sim.evalBondTension();
    glColor3f(1.0,0.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-1 );
    glColor3f(0.0,1.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.vel, 1e-1 );
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
}

SpaceCraftDynamicsOCLapp::SpaceCraftDynamicsOCLapp( int& id, int WIDTH_, int HEIGHT_ ) : SpaceCraftDynamicsApp( id, WIDTH_, HEIGHT_ ) { 

}

void SpaceCraftDynamicsOCLapp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    // Render simulation
    glLineWidth(0.5); 
    if(bDouble){
        SpaceCraftDynamicsApp::drawSim( *_sim    );
    }else{
        drawSim( *_sim_cl );
        //runSim_cpu( sim_cl );
    }

    //sim2.run_Cholesky(1);
    //sim2.run_LinSolve(1);
    //renderTruss( sim2.nBonds, sim2.bonds, sim2.points, sim2.strain, 1000.0 );

    // draw ring nodes
    // glLineWidth(3.0);
    // glColor3f(1.0,0.0,1.0);
    // for(const Ring* o: theSpaceCraft->rings){
    //     Draw3D::drawMatInPos( o->pose.rot, o->pose.pos, {100.,100.,100.} );
    //     Node** nds = (Node**)&o->nodes;
    //     for( int i=0; i<4; i++ ){
    //         //printf( "nds[%i] %li \n", i, (long)nds[i] );
    //         if( nds[i] == 0 ) continue;
    //         Draw3D::drawPointCross( nds[i]->pos, 10.0 );
    //     }
    // }

    // picker.hray = (Vec3d)(cam.rot.c);
    // picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    // glLineWidth(5.0);
    // if     (picker.edit_mode == EDIT_MODE::vertex){ if( picker.picked>=0 ){ Vec3f p = *(Vec3f*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawPointCross( p, 10.0 );                              } }
    // else if(picker.edit_mode == EDIT_MODE::edge  ){ if( picker.picked>=0 ){ Vec2i b = *(Vec2i*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawLine      ( sim_cl.points[b.x].f, sim_cl.points[b.y].f ); } }

};

// void SpaceCraftDynamicsOCLapp::drawHUD(){
//     glDisable( GL_LIGHTING );
//     glDisable(GL_DEPTH_TEST);
//     //gui.draw();
//     //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
//     //glPopMatrix();
// }

// ===================== MAIN

int main(int argc, char *argv[]){


    sim_cl.print_devices();
    sim_cl.initOCL();

    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	SpaceCraftDynamicsOCLapp * app = new SpaceCraftDynamicsOCLapp( junk, dm.w-150, dm.h-100 );
    app->_mesh   = &mesh;
    app->_sim    = &sim2;
    app->_sim_cl = &sim_cl;    

    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    printf( "argc %i \n", argc );
    LambdaDict funcs;
    funcs["-s"]={1,[&](const char** ss){ 
        reloadShip( ss[0], mesh );
        app->initSim( mesh );
    }}; 
    process_args( argc, argv, funcs );
	if( (sim2.nPoint==0) && (sim_cl.nPoint==0) ){ app->initSimDefault(); }
	app->loop( 1000000 );
	return 0;
}
















