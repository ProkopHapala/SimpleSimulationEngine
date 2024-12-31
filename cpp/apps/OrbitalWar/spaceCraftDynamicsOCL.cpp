
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

// ======================  Global Variables & Declarations

using namespace SpaceCrafting;

bool bRun = false;
bool bViewPointLabels = false;

Mesh::Builder2 mesh;
OrbSim         sim2;   // CPU double precision
OCL_Orb        sim_cl; // single precision both CPU and GPU

// make list of fixed points containing following indices : [0]
std::vector<int> fixPoints{ 1 };
double Kfix = 1e+12;

int glo_truss=0, glo_capsula=0, glo_ship=0;
double elementSize  = 5.;


// Render 
void runSim_gpu( OCL_Orb& sim, int niter=100 ){
    if( (sim.nPoint==0)||(sim.points==0) ){ printf( "ERROR in runSim() sim.nPoint=%i sim.points=%p => exit() \n", sim.nPoint, sim.points ); exit(0); };
    niter=1;
    //niter=10;
    int nbig = 20;
    int nsub = 20;
    niter = nbig*nsub;
    //sim.dt = 0.1;
    //sim.dt = 0.1e-4;

    //sim.set_time_step(0.1   );
    //sim.set_time_step(0.01  );
    //sim.set_time_step(0.001 );
    sim.set_time_step(0.0001);

    //niter = 500;
    if(bRun){
        long t0 = getCPUticks();
        //sim.run( niter, 1e-4, 1e-4 );
        //sim.run_omp( niter, false, 1e-3, 1e-4 );
        //sim.run_omp( niter, true, 1e-3, 1e-4 );
        //sim.damping = 1e-5; sim.dt      = 1e-3; // does not have any effect here
        //sim.run_ocl( niter*nsub, 0b001, 0b111 );   // 0b001 - upload points, 0b111 - download points, velocities, forces
        //sim.run_PDcl( nbig, nsub, 0b001, 0b111 ); // 0b001 - upload points, 0b111 - download points, velocities, forces

        sim.run_PDcl( nbig, nsub, 0b111, 0b111 );

        // ---- Debug loop --- check if forces are the same between GPU and CPU
        //printf( "runSim() fmax=%g fmax_ref=%g\n", sim.getFmax(), fmax_ref );
        //printf( "runSim() fmax=%g fmax_ref=%g\n", sim.getFmax(), fmax_ref );
        // for(int itr=0; itr<niter; itr++){
        //     sim.damping = 1e-5;
        //     sim.dt      = 1e-3;
        //     sim.run_ocl( 1, 0b011, 0b111 );
        //     //sim.cleanForce();
        //     //sim.evalTrussForces_neighs2();
        //     //sim.applyForceCentrifug( sim.rot0.f, sim.omega.f, sim.omega.w );
        //     sim.move_MD( 1e-3, 1e-5 );
        // }
        double T = (getCPUticks()-t0)*1e-6;
        //printf( "runSim() DONE T=%g[ms] %g[ms/iter] niter=%i,nP=%i,nE=%i \n", T, T/niter, niter, sim.nPoint, sim.nNeighMax );
        //printf( "runSim() nPoint=%i nBonds=%i \n", sim.nPoint, sim.nBonds );
    }
    sim.evalBondTension();
    glColor3f(1.0,0.0,1.0);
    //renderPoinSizes( sim.nPoint, sim.points, 0.001 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-6 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-3 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-2 );
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-1 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1.0 );

    glColor3f(0.0,1.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.vel, 1e-1 );
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
}

void runSim_cpu( OrbSim_f& sim, int niter=100 ){
    if( (sim.nPoint==0)||(sim.points==0) ){ printf( "ERROR in runSim_cpu() sim.nPoint=%i sim.points=%p => exit() \n", sim.nPoint, sim.points ); exit(0); };
    niter=1;
    //niter=10;
    int nbig = 1;
    int nsub = 20;
    niter = nbig*nsub;

    //niter = 500;
    if(bRun){
        long t0 = getCPUticks();
        sim.run_Cholesky_omp_simd(nbig);
        double T = (getCPUticks()-t0)*1e-6;
    }
    sim.evalBondTension();
    glColor3f(1.0,0.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-1 );
    glColor3f(0.0,1.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.vel, 1e-1 );
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
}

void runSim_double( OrbSim& sim, int niter=100 ){
    //printf( "runSim_double() nPoint=%i nBonds=%i nNeighMax=%i nNeighMaxLDLT=%i\n", sim.nPoint, sim.nBonds, sim.nNeighMax, sim.nNeighMaxLDLT );
    if( (sim.nPoint==0)||(sim.points==0) ){ printf( "ERROR in runSim_double() sim.nPoint=%i sim.points=%p => exit() \n", sim.nPoint, sim.points ); exit(0); };
    niter=1;
    //niter=10;
    int nbig = 1;
    int nsub = 20;
    niter = nbig*nsub;

    //niter = 500;
    if(bRun){
        long t0 = getCPUticks();
        sim.run_Cholesky_omp_simd(nbig);
        double T = (getCPUticks()-t0)*1e-6;
    }
    sim.evalBondTension();
    glColor3f(1.0,0.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-1 );
    glColor3f(0.0,1.0,1.0);
    renderPointForces( sim.nPoint, sim.points, sim.vel, 1e-1 );
    renderTruss      ( sim.nBonds, sim.bonds,  sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
}


// ======================  Free Functions

void to_OrbSim( OrbSim& sim, Mesh::Builder2& mesh, Vec3d p0, Vec3d ax ){
    //int nneighmax_min = 16;
    int nneighmax_min = 8;
    exportSim( sim, mesh, workshop,  nneighmax_min );
    for(int i=0; i<sim2.nPoint; i++)  sim2.points[i].w=1.0;
    for(int i=0; i<sim2.nBonds;  i++) sim2.params[i].y=10000.0;
    //sim2.printAllNeighs();

    sim.reallocFixed();
    for(int i:fixPoints){ printf("to_OrbSim fixing point %i @kFix=%p\n", i, sim.kFix); sim.kFix[i]=Kfix; }

    int n = sim.nPoint;

    double dt = 0.0001;
    sim.set_time_step( dt );

    mat2file<int>( "neighs_before.log",  n, sim.nNeighMax,      sim.neighs,     "%5i " );
    //sim2.prepare_Cholesky( 0.05, 32 );
    sim.prepare_LinearSystem( dt, true, true, true, 32 );
    mat2file<int>( "neighs_after.log",   n, sim.nNeighMaxLDLT,  sim.neighsLDLT, "%5i " );

    mat2file<double>( "PDmat.log",  n,n, sim.PDmat  );
    mat2file<double>( "LDLT_L.log", n,n, sim.LDLT_L );

    double omega = 0.0;
    sim.cleanVel();
    sim.addAngularVelocity(  p0, ax*omega );
    //apply_torq( sim.nPoint, p0, ax*omega, sim.points, sim.vel );  
}

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

void reloadShip( const char* fname  ){
    theSpaceCraft->clear();                  // clear all components
    //luaL_dofile(theLua, "data/spaceshil1.lua");
    printf("#### START reloadShip('%s')\n", fname );
    if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
    printf( "Lua::dofile(%s) DONE \n", fname );
    theSpaceCraft->checkIntegrity();

    mesh.clear();
    BuildCraft_truss( mesh, *theSpaceCraft, 30.0 );
    mesh.printSizes();

    to_OCL_Orb(sim_cl, mesh, Vec3fZero, Vec3fZ );
    to_OrbSim( sim2,   mesh, Vec3dZero, Vec3dZ );
    
    printf("#### END reloadShip('%s')\n", fname );
};

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

void makeTrussShape( int ishape, int nseg, double R, double r, bool bDouble=true, bool bSingle=true ){

    StickMaterial *o = new StickMaterial();
    //Material{ name="Kevlar", density=1.44e+3, Spull=3.6e+9, Spush=0.0,    Kpull=154.0e+9, Kpush=0.0,      reflectivity=0.6,  Tmelt=350 }
    //Material{ name="Steel" , density=7.89e+3, Spull=1.2e+9, Spush=1.2e+9, Kpull=200.0e+9, Kpush=200.0e+9, reflectivity=0.85, Tmelt=800 }
    //st1  = StickMaterial( "GS1_long", "Steel", 0.1,  0.005 )
    //st2  = StickMaterial( "GS1_perp", "Steel", 0.05, 0.003 )
    //st3  = StickMaterial( "GS1_in",   "Steel", 0.04, 0.002 )
    //st4  = StickMaterial( "GS1_out",  "Steel", 0.04, 0.002 )

    workshop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
    workshop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );

    Vec3d p0{0.0,0.0,0.0};
    Vec3d p1{R  ,0.0,0.0};
    Vec3d ax{0.0,0.0,1.0};
    Vec3d up{0.0,1.0,0.0};
    
    //BuildCraft_truss( mesh, *theSpaceCraft, 30.0 );
    mesh.clear();
    //mesh.block();
    //mesh.wheel( p0, p1, ax, nseg, 0.2 );
    //wheel( mesh, p0, p1, ax, nseg, Vec2d{0.2,0.2}, Quat4i{0,0,0,0} );
    switch(ishape){
        case 0: ngon   ( mesh, p0, p1, ax, nseg, 0 ); break;
        case 1: wheel  ( mesh, p0, p1, ax, nseg, Vec2d{r,r}, Quat4i{0,0,0,0}       ); break;
        case 2: girder1( mesh, p0, (p1-p0).normalized()*r*4*nseg, ax, nseg, r,          Quat4i{0,0,0,0}, true ); break;
        case 3: triangle_strip( mesh, p0, (p1-p0).normalized()*r*nseg, up, nseg, r, 0, true );
        //int girder1( Builder2& mesh, Vec3d p0, Vec3d p1, Vec3d up, int n, double width, Quat4i stickTypes ){
    }
    //wheel( mesh, o->pose.pos, o->pose.pos+o->pose.rot.b*o->R, o->pose.rot.c, o->nseg, o->wh, o->st );
    //Quat4i& b = mesh.blocks.back();
    //o->pointRange = {b.x,(int)mesh.verts.size()};
    //o->stickRange = {b.y,(int)mesh.edges.size()};
    mesh.printSizes();

    if(bDouble){ to_OrbSim ( sim2,   mesh,        p0,         ax      ); }
    if(bSingle){ to_OCL_Orb( sim_cl, mesh, (Vec3f)p0, (Vec3f)(ax*5.0) ); }

    distort_points( sim_cl.nPoint, sim_cl.points, sim_cl.points, 1.0, Vec3fOne, 15454 ); 

    // std::vector<Quat4f> ps_bak(sim_cl.nPoint); 
    // distort_points( sim_cl.nPoint, sim_cl.points, ps_bak.data(), 1.0, Vec3fOne, 15454 ); 

    // sim_cl.bmix.x = 1.0; test_SolverConvergence( mesh, ps_bak.data(), 3,   8, 0.5, 0.05 );
    // sim_cl.bmix.x = 1.0; test_SolverConvergence( mesh, ps_bak.data(), 5,   8, 0.5, 0.05 );
    // sim_cl.bmix.x = 1.0; test_SolverConvergence( mesh, ps_bak.data(), 10,  8, 0.5, 0.05 );
    // sim_cl.bmix.x = 1.0; test_SolverConvergence( mesh, ps_bak.data(), 20,  8, 0.5, 0.05 );
    // sim_cl.bmix.x = 1.0; test_SolverConvergence( mesh, ps_bak.data(), 50,  8, 0.5, 0.05 );
    // sim_cl.bmix.x = 1.0; test_SolverConvergence( mesh, ps_bak.data(), 100, 8, 0.5, 0.05 );

    printf("#### END makeShip_Whee()\n" );
    //exit(0);
};



// ====================== Class Definitions

class SpaceCraftEditorApp : public AppSDL2OGL_3D { public:

    PickerUI  picker;

    // ==== function declarations

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

	SpaceCraftEditorApp( int& id, int WIDTH_, int HEIGHT_ , int argc, char *argv[]);

};

SpaceCraftEditorApp::SpaceCraftEditorApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //Lua1.init();
    fontTex       = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    Draw::fontTex = fontTex;

    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
    if(argc<=1){
        //reloadShip( "data/ship_ICF_interceptor_1.lua" );
        //makeTrussShape( 2, 1,   100.0, 10.0,   true, true );
        //makeTrussShape( 2, 100, 100.0, 10.0,  true, true );
        makeTrussShape  ( 2, 10, 100.0, 10.0,  true, true );
        //makeTrussShape  ( 3, 50, 100.0, 10.0,  true, true );
    }

    picker.picker = &sim_cl;   picker.Rmax=10.0;

    VIEW_DEPTH = 10000.0;
    zoom       = 1000.0;
    printf( "### SpaceCraftEditorApp() DONE\n" );
}

void SpaceCraftEditorApp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    // Render simulation
    glLineWidth(0.5); 
    bool bDouble = true;
    if(bDouble){
        runSim_double( sim2  );
    }else{
        runSim_gpu( sim_cl );
        //runSim_cpu( sim_cl );
    }

    //sim2.run_Cholesky(1);
    //sim2.run_LinSolve(1);
    //renderTruss( sim2.nBonds, sim2.bonds, sim2.points, sim2.strain, 1000.0 );

    // draw ring nodes
    glLineWidth(3.0);
    glColor3f(1.0,0.0,1.0);
    for(const Ring* o: theSpaceCraft->rings){
        Draw3D::drawMatInPos( o->pose.rot, o->pose.pos, {100.,100.,100.} );
        Node** nds = (Node**)&o->nodes;
        for( int i=0; i<4; i++ ){
            //printf( "nds[%i] %li \n", i, (long)nds[i] );
            if( nds[i] == 0 ) continue;
            Draw3D::drawPointCross( nds[i]->pos, 10.0 );
        }
    }

    picker.hray = (Vec3d)(cam.rot.c);
    picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    glLineWidth(5.0);
    if     (picker.edit_mode == EDIT_MODE::vertex){ if( picker.picked>=0 ){ Vec3f p = *(Vec3f*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawPointCross( p, 10.0 );                              } }
    else if(picker.edit_mode == EDIT_MODE::edge  ){ if( picker.picked>=0 ){ Vec2i b = *(Vec2i*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawLine      ( sim_cl.points[b.x].f, sim_cl.points[b.y].f ); } }

};

void SpaceCraftEditorApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    //gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}


void SpaceCraftEditorApp::keyStateHandling( const Uint8 *keys ){
    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
    if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.roll(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.roll( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }

    if( keys[ SDL_SCANCODE_LEFTBRACKET  ] ){ theSpaceCraft->nodes[7]->calong-=0.001; theSpaceCraft->nodes[7]->updateBound(); }
    if( keys[ SDL_SCANCODE_RIGHTBRACKET ] ){ theSpaceCraft->nodes[7]->calong+=0.001; theSpaceCraft->nodes[7]->updateBound(); }
};

void SpaceCraftEditorApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    //if(event.type == SDL_MOUSEWHEEL){
    //    if     (event.wheel.y > 0){ zoom*=VIEW_ZOOM_STEP; }
    //    else if(event.wheel.y < 0){ zoom/=VIEW_ZOOM_STEP; }
    //}
    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    //gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_m: picker.switch_mode(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_l:
                //    //reloadShip( );
                //    //onSelectLuaShipScript.GUIcallback(lstLuaFiles);
                //    break;
                case SDLK_SPACE: bRun = !bRun; break;
                case SDLK_l:     bViewPointLabels ^=1; break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  picker.pick(); break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT: break;
                case SDL_BUTTON_RIGHT:break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

    //printf( "compGui %li \n", compGui );
    //if(compGui) if( compGui->check() ){ renderShip(); }

}

// ===================== MAIN

LambdaDict funcs;
SpaceCraftEditorApp * app;

int main(int argc, char *argv[]){

    printf( "argc %i \n", argc );
    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    //funcs["-s"]={1,[&](const char** ss){ app->reloadShip( ss[0] ); }}; 
    funcs["-s"]={1,[&](const char** ss){ reloadShip( ss[0] ); }}; 

    sim_cl.print_devices();
    sim_cl.initOCL();
    //sim_cl.initOCL(0);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	// https://www.opengl.org/discussion_boards/showthread.php/163904-MultiSampling-in-SDL
	//https://wiki.libsdl.org/SDL_GLattr
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
    glEnable(GL_MULTISAMPLE);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	app = new SpaceCraftEditorApp( junk , dm.w-150, dm.h-100, argc, argv );
    process_args( argc, argv, funcs );
	//app = new SpaceCraftEditorApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















