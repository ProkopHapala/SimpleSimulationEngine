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
#include "TrussDynamics_d.h"
#include "TrussDynamics_f.h"
//#include "OCL_Orb.h"
#include "spaceCraftSimulator.h"

#include "testUtils.h"

namespace SpaceCrafting {

// =======================================================================
class SpaceCraftDynamicsApp : public AppSDL2OGL_3D { public:

    SpaceCraftSimulator* simulator=0;
    Mesh::Builder2* _mesh=0;
    TrussDynamics_d*       _sim =0;
    TrussDynamics_f*       _sim_f=0;

    double Estrain = 0.0;

    // === View Options
    bool bDouble             = true;
    bool bRun                = false;
    bool bViewPointLabels    = false;
    bool bViewFixedPoints    = true;
    bool bDrawTrj            = false;
    bool bViewStrainNumbers  = false;
    bool bViewMassNumbers    = false;
    bool bViewBondStiffness  = false;
    bool bViewResudualForces = false;
    bool bViewVelocities     = false;

    double scale_force       = 0.01;
    double scale_velocity    = 1.0;
    float fontSize3D         = 0.014;
    int perFrame = 1;
    //int perFrame = 10;
    //int perFrame = 100;

    // === GUI
    GUI gui;
    CheckBoxList* panel_view=0;
    MultiPanel*   panel_sim=0;

    // === Functions

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    virtual void initGUI();
    //virtual void mouseHandling( );

    virtual void bindSimulators( SpaceCraftSimulator* simulator_ ){
        simulator=simulator_;
        _mesh  = simulator->getMesh();
        _sim   = simulator->getTrussSim();
        _sim_f = simulator->getTrussSim_f();
    }

    void drawSim  ( TrussDynamics_d&   sim );
    void drawSim_f( TrussDynamics_f& sim );
    //virtual void initSimDefault();

	SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_ );

    

};

void SpaceCraftDynamicsApp::initGUI(){
    // === View Options Panel
    panel_view = new CheckBoxList(5, 5, 200, fontSizeDef*2);
    gui.addPanel(panel_view);
    panel_view->addBox("Point Labels",     &bViewPointLabels);
    panel_view->addBox("Fixed Points",     &bViewFixedPoints);
    panel_view->addBox("Trajectories",     &bDrawTrj);
    panel_view->addBox("Strain Numbers",   &bViewStrainNumbers);
    panel_view->addBox("Mass Numbers",     &bViewMassNumbers);
    panel_view->addBox("Bond Stiffness",   &bViewBondStiffness);
    panel_view->addBox("Residual Forces",  &bViewResudualForces);
    panel_view->addBox("Velocities",       &bViewVelocities);

    // === Simulation Parameters Panel
    panel_sim = new MultiPanel("Simulation", 210, 5, 400, fontSizeDef*2, 4);
    gui.addPanel(panel_sim);
    panel_sim->subs[0]->setValue(scale_force)->setRange(0.001, 0.1);
    panel_sim->subs[1]->setValue(scale_velocity)->setRange(0.1, 10.0);
    panel_sim->subs[2]->setValue(fontSize3D)->setRange(0.001, 0.1);
    panel_sim->subs[3]->setValue(perFrame)->setRange(1, 100);
}

void SpaceCraftDynamicsApp::drawSim( TrussDynamics_d& sim ){
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    glColor3f( 0.0f,0.0f, 0.0f );
    //if(bViewResudualForces){ glColor3f( 1.0f,0.0f, 0.0f ); Draw3D::drawVectorArray( sim.nPoint, sim.ps_cor, sim.linsolve_b, scale_force    );  }
    if(bViewMassNumbers    ){ drawPointProperties( sim.nPoint, sim.points, (double*)sim.points, 4,3, fontTex,  fontSize3D, "%.0f" ); };
    if(bViewBondStiffness  ){ drawBondProperties ( sim.nBonds, sim.bonds, sim.points, (double*)sim.bparams, 4,2, fontTex,  fontSize3D, "%.3e"); };
    if(bViewStrainNumbers  ){ drawBondProperties ( sim.nBonds, sim.bonds, sim.points, sim.strain, 1,0, fontTex,  fontSize3D, "%.3f", 100.0 ); };
    if(bViewResudualForces ){ glColor3f( 1.0f,0.0f, 0.0f ); renderPointForces      ( sim.nPoint, sim.points, sim.forces,     scale_force );     }
    if(bViewVelocities     ){ glColor3f( 0.0f,0.5f, 0.0f ); renderPointForces      ( sim.nPoint, sim.points, sim.vel,        scale_velocity );  }
    if(bViewPointLabels    ){ pointLabels( sim.nPoint, sim.points, fontTex, fontSize3D );}
    if(bViewFixedPoints && (sim.kFix!=0) ){ renderPoinsSizeRange( sim.nPoint, sim.points, sim.kFix, Vec2d{ 1.0, 1e+300 }, 10.0 ); }
};

void SpaceCraftDynamicsApp::drawSim_f( TrussDynamics_f& sim ){
    Estrain = sim.evalBondTension();
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    glColor3f( 0.0f,0.0f, 0.0f );
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
    initGUI();
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
            Estrain = _sim->evalBondTension();
        }
        drawSim( *_sim   );
    }else{
        if(bRun){
            //_sim_f->run_Cholesky_omp_simd(perFrame);
            _sim_f->run_LinSolve( perFrame );
            Estrain = _sim_f->evalBondTension();
        }
        drawSim_f( *_sim_f );
    }
    double T = (getCPUticks()-t0);
    if(bRun)printf( "SpaceCraftDynamicsApp::drawSim(bDouble=%i,method=%i,bmix=%g,dt=%g,nsolve=%i) perFrame: %3i nPoint:%6i TIME: %8.3f [Mticks] %8.1f [tick/point]  Estrain: %10.2e\n", 
                                     bDouble, _sim->linSolveMethod, _sim->mixer.b_end, _sim->dt, _sim->nSolverIters, perFrame, _sim->nPoint, T*1e-6,  T/(perFrame*_sim->nPoint), Estrain );
	//if(!bDrawTrj)glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	//glEnable(GL_DEPTH_TEST);
    Draw3D::drawAxis( 1.0 );
};

void SpaceCraftDynamicsApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);
    gui.draw();
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
    if( gui.onEvent( mouseX, mouseY, event ) )return;

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
