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
#include "OrbSim_f.h"
//#include "OCL_Orb.h"
#include "spaceCraftSimulator.h"

namespace SpaceCrafting {

// void reloadShip( const char* fname, Mesh::Builder2& mesh ){
//     printf("#### START reloadShip('%s')\n", fname );
//     theSpaceCraft->clear();                       DEBUG
//     //luaL_dofile(theLua, "data/spaceshil1.lua"); DEBUG
//     if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
//     printf( "Lua::dofile(%s) DONE \n", fname );
//     theSpaceCraft->checkIntegrity();
//     mesh.clear();
//     BuildCraft_truss( mesh, *theSpaceCraft, 30.0 );
//     mesh.printSizes();
//     printf("#### END reloadShip() \n");
// };

// void to_OrbSim( OrbSim& sim, Mesh::Builder2& mesh, double dt=0.1, Vec3d p0=Vec3dZero, Vec3d omega=Vec3dZero ){
//     //printf("#### START to_OrbSim() \n");
//     exportSim( sim, mesh, workshop );
//     std::vector<int> fixPoints{ 2 }; sim.setFixPoints( fixPoints.size(), fixPoints.data() );
//     sim.prepare_LinearSystem( dt, true, true, true, 256 );
//     //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CG;
//     sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Cholesky;
//     //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CholeskySparse;
//     sim.cleanVel();
//     if( omega.norm2()>1e-16 )sim.addAngularVelocity( p0, omega );
//     sim.addAngularVelocity( p0, Vec3dZ*0.01 );
// };

// void to_OrbSim_f( OrbSim_f& sim, Mesh::Builder2& mesh, double dt=0.1, Vec3f p0=Vec3fZero, Vec3f omega=Vec3fZero ){
//     //printf("#### START to_OrbSim() \n");
//     exportSim( sim, mesh, workshop );
//     if( sim.nPoint==0 ){ printf( "ERROR in to_OrbSim() sim.nPoint=%i => exit() \n", sim.nPoint ); exit(0); };
//     //std::vector<int> fixPoints{ 2 }; sim.setFixPoints( fixPoints.size(), fixPoints.data() );
//     sim.prepare_LinearSystem( dt, true, true, true, 256 );
//     //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CG;
//     sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Cholesky;
//     //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::CholeskySparse;
//     sim.cleanVel();
//     if( omega.norm2()>1e-16 )sim.addAngularVelocity( p0, omega );
//     sim.addAngularVelocity( p0, Vec3fZ*0.01 );
// };

// void init_workshop(){
//     //StickMaterial *o = new StickMaterial();
//     workshop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
//     workshop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );
// }

// =======================================================================
class SpaceCraftDynamicsApp : public AppSDL2OGL_3D { public:

    SpaceCraftSimulator* simulator=0;
    Mesh::Builder2* _mesh=0;
    OrbSim*         _sim =0;
    OrbSim_f*       _sim_f=0;

    bool bDouble          = false;
    bool bRun             = false;
    bool bViewPointLabels = false;
    bool bViewFixedPoints = true;
    bool bDrawTrj         = false;

    int perFrame = 1;
    //int perFrame = 10;
    //int perFrame = 100;
    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

    virtual void bindSimulators( SpaceCraftSimulator* simulator_ ){
        simulator=simulator_;
        _mesh  = simulator->getMesh();
        _sim   = simulator->getOrbSim();
        _sim_f = simulator->getOrbSim_f();
    }

    void drawSim  ( OrbSim&   sim );
    void drawSim_f( OrbSim_f& sim );
    //virtual void initSimDefault();

	SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_ );

};

// void SpaceCraftDynamicsApp::initSimDefault( ){
//     init_workshop ();
//     //makeTrussShape( *_mesh, 3, 10, 100.0, 10.0 );
//     makeTrussShape( *_mesh, 2, 100, 100.0, 10.0 );
//     to_OrbSim  ( *_sim,   *_mesh );
//     to_OrbSim_f( *_sim_f, *_mesh );
// }

void SpaceCraftDynamicsApp::drawSim( OrbSim& sim ){
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
    if(bViewFixedPoints && (sim.kFix!=0) ) renderPoinsSizeRange( sim.nPoint, sim.points, sim.kFix, Vec2d{ 1.0, 1e+300 }, 10.0 );
};

void SpaceCraftDynamicsApp::drawSim_f( OrbSim_f& sim ){
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    if(bViewPointLabels) pointLabels( sim.nPoint, sim.points, fontTex, 0.02 );
    if(bViewFixedPoints && (sim.kFix!=0) ) renderPoinsSizeRange( sim.nPoint, sim.points, sim.kFix, Vec2f{ 1.0, 1e+300 }, 10.0f );
};

SpaceCraftDynamicsApp::SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //Lua1.init();
    fontTex     = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
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
    
    long t0 = getCPUticks(); 
    if( bDouble ){
        if(bRun){
            //_sim->run_Cholesky_omp_simd(perFrame);
            _sim->run_LinSolve( perFrame );
        }
        drawSim( *_sim   );
    }else{
        if(bRun){
            //_sim_f->run_Cholesky_omp_simd(perFrame);
        }
        drawSim_f( *_sim_f );
    }
    double T = (getCPUticks()-t0);
    if(bRun)printf( "SpaceCraftDynamicsApp::drawSim(bDouble=%i) perFrame: %3i nPoint:%6i TIME: %8.3f [Mticks] %8.1f [tick/point] \n", bDouble, perFrame, _sim->nPoint, T*1e-6,  T/(perFrame*_sim->nPoint) );
	//if(!bDrawTrj)glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	//glEnable(GL_DEPTH_TEST);
    Draw3D::drawAxis( 1.0 );
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
