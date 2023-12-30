
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

#include "OrbSim.h"
#include "OCL_Orb.h"
//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

//#include "Tree.h"
//#include "spaceCraftEditorUtils.h"

#include "SpaceCraftGUI.h"
#include "argparse.h"

// ======================  Global Variables & Declarations

using namespace SpaceCrafting;

bool bRun = false;

Mesh::Builder2 mesh2;
//OrbSim_f       sim;
OCL_Orb        sim;

int glo_truss=0, glo_capsula=0, glo_ship=0;
double elementSize  = 5.;

// Render 
void runSim( OCL_Orb& sim, int niter=100 ){
    niter=1;
    //niter=10;
    niter = 100;
    //niter = 500;
    if(bRun){
        long t0 = getCPUticks();
        //sim.run( niter, 1e-4, 1e-4 );
        //sim.run_omp( niter, false, 1e-3, 1e-4 );
        //sim.run_omp( niter, true, 1e-3, 1e-4 );
        //sim.damping = 1e-5; sim.dt      = 1e-3; // does not have any effect here
        sim.run_ocl( niter, 0b001, 0b111 ); // 0b001 - upload points, 0b111 - download points, velocities, forces
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
        printf( "runSim() DONE T=%g[ms] %g[ms/iter] niter=%i,nP=%i,nE=%i \n", T, T/niter, niter, sim.nPoint, sim.nNeighMax );
        //printf( "runSim() nPoint=%i nBonds=%i \n", sim.nPoint, sim.nBonds );
    }
    sim.evalBondTension();
    //renderPoinSizes( sim.nPoint, sim.points, 0.001 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-6 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-3 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-2 );
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-1 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1.0 );
}

// ======================  Free Functions

void reloadShip( const char* fname  ){
    theSpaceCraft->clear();                  // clear all components
    //luaL_dofile(theLua, "data/spaceshil1.lua");
    printf("#### START reloadShip('%s')\n", fname );
    Lua::dofile(theLua,fname);
    printf( "Lua::dofile(%s) DONE \n", fname );
    theSpaceCraft->checkIntegrity();

    mesh2.clear();
    BuildCraft_truss( mesh2, *theSpaceCraft, 30.0 );
    mesh2.printSizes();
    exportSim( sim, mesh2, workshop );

    // OpenCL initialization
    printf("###### OpenCL initialization\n");
    sim.makeKrenels_Orb( "./common_resources/cl" );
    sim.initCLBuffsOrb(  );
    
    //sim.setup_test_enque();
    //sim.test_enque();

    sim.damping = 1e-5; sim.dt = 1e-3; // must be here before sim.setup_evalTrussForce2();
    sim.setup_evalTrussForce1();
    sim.setup_evalTrussForce2();
    sim.setup_move();
    sim.setup_assembleAndMove();
    sim.setup_evalTrussBondForce();
    
    printf("#### END reloadShip('%s')\n", fname );
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
    if(argc<=1)reloadShip( "data/ship_ICF_interceptor_1.lua" );

    picker.picker = &sim;   picker.Rmax=10.0;

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

    runSim( sim );
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
 
    picker.hray = (Vec3d)(cam.rot.c);
    picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    glLineWidth(5.0);
    if     (picker.edit_mode == EDIT_MODE::vertex){ if( picker.picked>=0 ){ Vec3f p = *(Vec3f*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawPointCross( p, 10.0 );                              } }
    else if(picker.edit_mode == EDIT_MODE::edge  ){ if( picker.picked>=0 ){ Vec2i b = *(Vec2i*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawLine      ( sim.points[b.x].f, sim.points[b.y].f ); } }

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
                case SDLK_l:
                    //reloadShip( );
                    //onSelectLuaShipScript.GUIcallback(lstLuaFiles);
                    break;
                case SDLK_SPACE: bRun = !bRun; break;
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

    sim.print_devices();
    sim.initOCL();
    //sim.initOCL(0);

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
















