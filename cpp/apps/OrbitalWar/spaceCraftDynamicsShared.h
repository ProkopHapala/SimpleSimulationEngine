#ifndef SpaceCraftDynamicsShared_h
#define SpaceCraftDynamicsShared_h

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
//#include "OCL_Orb.h"

namespace SpaceCrafting {

void reloadShip( const char* fname, Mesh::Builder2& mesh ){
    printf("#### START reloadShip('%s')\n", fname );
    theSpaceCraft->clear();                       DEBUG
    //luaL_dofile(theLua, "data/spaceshil1.lua"); DEBUG
    if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
    printf( "Lua::dofile(%s) DONE \n", fname );
    theSpaceCraft->checkIntegrity();
    mesh.clear();
    BuildCraft_truss( mesh, *theSpaceCraft, 30.0 );
    mesh.printSizes();
    printf("#### END reloadShip() \n");
};

void to_OrbSim( OrbSim& sim, Mesh::Builder2& mesh, double dt=0.02, Vec3d p0=Vec3dZero, Vec3d omega=Vec3dZero ){
    //printf("#### START to_OrbSim() \n");
    exportSim( sim, mesh, workshop );
    std::vector<int> fixPoints{ 2 }; sim.setFixPoints( fixPoints.size(), fixPoints.data() );
    sim.prepare_LinearSystem( dt, true, true, true, 256 );
    //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CG;
    sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Cholesky;
    //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CholeskySparse;
    sim.cleanVel();
    if( omega.norm2()>1e-16 )sim.addAngularVelocity( p0, omega );
};

void init_workshop(){
    //StickMaterial *o = new StickMaterial();
    workshop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
    workshop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );
}

// =======================================================================
class SpaceCraftDynamicsApp : public AppSDL2OGL_3D { public:

    Mesh::Builder2* _mesh=0;
    OrbSim*         _sim =0;

    bool bRun             = false;
    bool bViewPointLabels = false;
    bool bDrawTrj=false;
    double time=0;
    //int perFrame = 1;
    //int perFrame = 10;
    int perFrame = 100;
    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

    void drawSim( OrbSim& sim );

	SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] );

};

void SpaceCraftDynamicsApp::drawSim( OrbSim& sim ){
    //timeit( "TIME sim.run_LinSolve(perFrame) T = %g [ms]\n", 1e-6, [&](){sim.run_LinSolve(perFrame);});
    if(bRun){
    long t0 = getCPUticks();  sim.time_LinSolver=0; sim.time_cg_dot=0; sim.cgSolver_niterdone=0;
    //sim.run_LinSolve(perFrame);
    sim.run_Cholesky_omp_simd(perFrame);
    //sim.run_Cholesky_omp(perFrame);
    double time_run_LinSolve = (getCPUticks()-t0)*1e-6;
    printf( "TIME run: %g [Mticks] LinSolve: %g(%g\%) Mdot: %g(%g\%) niter_av=%g err2tot=%g \n", time_run_LinSolve, sim.time_LinSolver,   100*sim.time_LinSolver/time_run_LinSolve, sim.time_cg_dot, 100*sim.time_cg_dot/time_run_LinSolve, sim.cgSolver_niterdone/((double)perFrame), sim.cgSolver.err2tot );
    }
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
};

SpaceCraftDynamicsApp::SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //Lua1.init();
    fontTex     = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
    if(argc<=1){
        init_workshop ();
        makeTrussShape( *_mesh, 3, 10, 100.0, 10.0 );
        to_OrbSim( *_sim, *_mesh );
    }
    //exit(0);
    VIEW_DEPTH = 10000.0;
    zoom = 10.0;
}

void SpaceCraftDynamicsApp::draw(){
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
    drawSim( *_sim );
	//if(!bDrawTrj)glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	//glEnable(GL_DEPTH_TEST);
};

void SpaceCraftDynamicsApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);
    //gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

void SpaceCraftDynamicsApp::keyStateHandling( const Uint8 *keys ){
    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }
};

void SpaceCraftDynamicsApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    //if(event.type == SDL_MOUSEWHEEL){
    //    if     (event.wheel.y > 0){ zoom*=VIEW_ZOOM_STEP; }
    //    else if(event.wheel.y < 0){ zoom/=VIEW_ZOOM_STEP; }
    //}

    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_t: bDrawTrj=!bDrawTrj; glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); break;
                case SDLK_m:  break;
                //case SDLK_h:  warrior1->tryJump(); break;
                case SDLK_SPACE: bRun             ^=1; break;
                case SDLK_l:     bViewPointLabels ^=1; break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

}

}

#endif // SpaceCraftDynamicsShared_h
