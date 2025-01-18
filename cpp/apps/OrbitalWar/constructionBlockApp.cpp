
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "globals.h"

// Math & Geometry
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Solids.h"
#include "raytrace.h"

// Graphics & GUI
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"

#include "ConstructionBlock.h"
#include "MeshBuilder2.h"
#include "MeshBuilder2Draw.h"
#include "argparse.h"


//#include "testUtils.h"

// ======================  Global Variables & Declarations
Mesh::Builder2 mesh2;
ConstructionBlock block;

// ====================== Class Definitions

class ConstructionBlockApp : public AppSDL2OGL_3D { public:

    int perFrame = 10;

	bool bSmoothLines = 1;
	bool bWireframe   = 1;

	//DropDownList lstLuaFiles;
    GUI gui;

    //EDIT_MODE edit_mode = EDIT_MODE::component;
    //int picked = -1;
    //Vec3d hray; //= (Vec3d)(cam.rot.c);
    //Vec3d ray0; //= (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);

    int picked_block = -1;
    int ipick = -1;

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

	ConstructionBlockApp( int& id, int WIDTH_, int HEIGHT_ , int argc, char *argv[]);

};

ConstructionBlockApp::ConstructionBlockApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //Lua1.init();
    fontTex       = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    Draw::fontTex = fontTex;

    //block.faces[0].typ=1;
    block.Ls=Vec3d{1.1,1.0,0.9};
    for(int i=0;i<6;i++){
        //block.faces[i].typ=2;
        block.faces[i].typ=3;
    }
    drawBlock( mesh2, block );

    //printf( "mesh2.tris.size(): \n", mesh2.tris.size() );

    //plateGui  = (PlateGUI* )gui.addPanel( new PlateGUI ( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );
    //girderGui = (GirderGUI*)gui.addPanel( new GirderGUI( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );

    //truss.loadXYZ(  "data/octShip.xyz" );
    //DropDownList* lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    //lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    //lstLuaFiles->setCommand([this](GUIAbstractPanel* panel){ onSelectLuaShipScript.GUIcallback(panel); });

}


void ConstructionBlockApp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
    glLightfv( GL_LIGHT0, GL_POSITION,  (float*)&cam.rot.c  );

    glColor3f( 1.0,1.0,1.0 );
    drawFaces( mesh2 );

    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    // Draw All Edges
    glColor3f(0.0,0.0,0.0);
    glLineWidth(1.0);
    drawEdges( mesh2 );

    if(ipick>=0){
        //printf( "ipick %i \n", ipick );
        if( mesh2.selection_mode == (int)Mesh::Builder2::SelectionMode::face ){
            Quat4i ch = mesh2.chunks[ipick];
            glLineWidth(5.0);
            glColor3f(0.0,0.7,0.0);
            glBegin(GL_LINE_LOOP);
            for(int i=0;i<ch.z;i++){
                int iv = mesh2.strips[ch.x+i];
                //printf( " %i ", i, ch.x+i, iv );
                Draw3D::vertex( mesh2.verts[iv].pos );
            }
            glEnd();
        };
    }

    // Draw selected edges
    if( mesh2.selection_mode==(int)Mesh::Builder2::SelectionMode::edge ){
        glColor3f(0.0,0.7,0.0);
        glLineWidth(5.0);
        drawSelectedEdges( mesh2 );
        drawSelectedEdgeLabels( mesh2, 0.02 );
    }

    glColor3f(0.f,0.f,0.f);
    drawPointLabels       ( mesh2, 0.02 );
    //drawEdgeLabels        ( mesh2, 0.02 );
    //drawSelectedEdgeLabels( mesh2, 0.02 );

    glLineWidth(5.0);
    Draw3D::drawAxis(10.0);

    glLineWidth(1.0);
    glColor3f(0.0,0.7,0.0);
    if(bDragging){ drawMuseSelectionBox(); }
};

void ConstructionBlockApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    // void Draw::drawText( const char * str, int itex, float sz, Vec2i block_size ){

    // sprintf(str_tmp, "time=%10.5f[s] mass=%g cog(%g,%g,%g) vcog(%g,%g,%g) L(%g,%g,%g) torq(%g,%g,%g) |F|=%g \n", sim.time, sim.mass, sim.cog.x,sim.cog.y,sim.cog.z, sim.vcog.x,sim.vcog.y,sim.vcog.z, sim.L.x,sim.L.y,sim.L.z, sim.torq.x,sim.torq.y,sim.torq.z, sim.F_residual );
    // //sprintf( str_tmp, "time= %10.5f[s] \n ", sim.time );
    // Draw::drawText( str_tmp, fontTex, fontSizeDef,  {WIDTH,HEIGHT-20}  );

    gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

void ConstructionBlockApp::keyStateHandling( const Uint8 *keys ){
    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3d){0.0,0.0,0.0} );
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
    //if( keys[ SDL_SCANCODE_LEFTBRACKET  ] ){ theSpaceCraft->nodes[7]->calong-=0.001; theSpaceCraft->nodes[7]->updateBound(); }
    //if( keys[ SDL_SCANCODE_RIGHTBRACKET ] ){ theSpaceCraft->nodes[7]->calong+=0.001; theSpaceCraft->nodes[7]->updateBound(); }



};

void ConstructionBlockApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){

         // ========  Mouse events  ========

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:{
                    //printf("SDL_BUTTON_LEFT: mouse_ray0 %g %g %g mouse_begin %g %g \n", ray0.x, ray0.y, ray0.z, mouse_begin_x, mouse_begin_y); 
                    mouseStartSelectionBox(); 
                    
                    //picker.pick();
                } break;

                case SDL_BUTTON_RIGHT:{} break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if( ray0.dist2(ray0_start)<0.1 ){ // too small for selection box 
                        ipick = mesh2.pickSelect( (Vec3d)ray0, (Vec3d)cam.rot.c, 0.1 );
                        printf( "ipick %i \n", ipick );
                    }else{
                        //ipick=-1;
                        mesh2.selectRect( (Vec3d)ray0_start, (Vec3d)ray0, (Mat3d)cam.rot );
                    }
                    bDragging=false;
                    break;
                case SDL_BUTTON_RIGHT: { 
                    ipick=-1;
                } break;
            }
            break;

        // ========  Keyboard events ========

        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_l:{} break;
                case SDLK_f:{
                    mesh2.selectionToFace();
                } break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

}

// ===================== MAIN

int main(int argc, char *argv[]){    
    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	ConstructionBlockApp * app = new ConstructionBlockApp( junk , dm.w-150, dm.h-100, argc, argv );
    //app->bindSimulators( &W ); 

    LambdaDict funcs;
    //funcs["-s"]={1,[&](const char** ss){ reloadShip( ss[0] ); }}; 
    //funcs["-perframe"]={1,[&](const char** ss){            sscanf( ss[0], "%i", &app->perFrame );              printf( "COMMAND LINE: -perframe(%i) \n", app->perFrame ); } };
    //funcs["-method"  ]={1,[&](const char** ss){ int im;    sscanf( ss[0], "%i", &im );  sim.linSolveMethod=im; printf( "COMMAND LINE: -method(%i)   \n", im            ); } };
    //funcs["-dt"      ]={1,[&](const char** ss){ float dt;  sscanf( ss[0], "%f", &dt );  sim.dt=dt;             printf( "COMMAND LINE: -dt( dt: %f ) \n", sim.dt        ); } };
    //funcs["-bmix"    ]={1,[&](const char** ss){ int istart; float bmix;  sscanf( ss[0], "%i,%f", &istart, &bmix ); W.sim.mixer.b_end=bmix; W.sim.mixer.istart=istart; printf( "COMMAND LINE: -bmix( istart:%i bmix: %f ) \n", W.sim.mixer.istart, W.sim.mixer.b_end );    } };
    //funcs["-fix"     ]={1,[&](const char** ss){ int n =  readlist( ss[0], W.fixPoints); printf("COMMAND LINE: -fix[%i]{%s}\n", n, ss[0] );  } };
    //funcs["-nsolve"  ]={1,[&](const char** ss){ int nsolv; sscanf( ss[0], "%i", &nsolv ); printf( "COMMAND LINE: -nsolve(%i) \n", nsolv ); W.sim_f.nSolverIters=nsolv; W.sim.nSolverIters=nsolv;  } };    

    process_args( argc, argv, funcs );

	//app = new ConstructionBlockApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
