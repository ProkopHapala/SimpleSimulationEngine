
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

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Truss.h"
#include "SpaceCraft.h"
#include "SpaceCraftDraw.h"
#include "SoftBody.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

#include "EditSpaceCraft.h"

#include "Tree.h"

#include "spaceCraftEditorUtils.h"

using namespace SpaceCrafting;

enum class EDIT_MODE : int { vertex=0, edge=1, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

//SpaceCraft craft;
Truss      truss;
int glo_truss=0, glo_capsula=0, glo_ship=0;

//char str[8096];

void reloadShip( const char* fname  ){
    theSpaceCraft->clear();
    //luaL_dofile(theLua, "data/spaceshil1.lua");
    printf("reloadShip('%s')", fname );
    Lua::dofile(theLua,fname);
    //luaL_dostring(theLua, "print('LuaDEBUG 1'); n1 = Node( {-100.0,0.0,0.0} ); print('LuaDEBUG 2'); print(n1)");

    theSpaceCraft->toTruss();

    if(glo_ship){ glDeleteLists(glo_ship,1); };
    glo_ship = glGenLists(1);
    glNewList( glo_ship, GL_COMPILE );
    drawSpaceCraft( *theSpaceCraft, 1 );
    glEndList();
};

class SpaceCraftEditGUI : public AppSDL2OGL_3D { public:

	class OnSelectLuaShipScript : public GUIEventCallback{ public:
        virtual int GUIcallback(GUIAbstractPanel* caller) override {
            ((DropDownList*)caller)->selectedToStr(str+sprintf(str,"data/"));
            reloadShip(str);
        }
    };

	bool bSmoothLines = 1;
	bool bWireframe   = 1;

	//DropDownList lstLuaFiles;
    GUI gui;
    DropDownList* lstLuaFiles;
    OnSelectLuaShipScript onSelectLuaShipScript;

    EDIT_MODE edit_mode = EDIT_MODE::vertex;
    //EDIT_MODE edit_mode = EDIT_MODE::edge;
    int picked = -1;
    Vec3d mouse_ray0;

    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

	SpaceCraftEditGUI( int& id, int WIDTH_, int HEIGHT_ );




};

SpaceCraftEditGUI::SpaceCraftEditGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //truss.loadXYZ(  "data/octShip.xyz" );
    //DropDownList* lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    //lstLuaFiles->onSelect = new OnSelectLuaShipScript();
    lstLuaFiles->onSelect = &onSelectLuaShipScript;

    TreeView* tvDir      = new TreeView    ( "DirView",20,HEIGHT_-400,200,20 ); gui.addPanel(tvDir);
    dir2tree(tvDir->root, "data" );
    tvDir->updateLines();

    /*
    char str[1000];
    getcwd( str, 1000);
    printf(" getcwd '%s'\n", str );
    */

    //std::vector<std::string> luaFiles;
    listDirContaining( "data", ".lua", lstLuaFiles->labels );

    //Lua1.init();
    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    Vec3f pf;
    Vec3d pd;
    //pd.set(1.5d,1.800001d,1.9d);
    pd.set(1.5d,1.800000001d,1.9d);
    //pd = pf;
    //pd = (Vec3d)pf;
    pf = (Vec3f)pd;
    print(pf);printf(" // pf\n");
    print(pd);printf(" // pd\n");
    //exit(0);

    //camera();

    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
    //reloadShip( "data/spaceshil1.lua" );
    onSelectLuaShipScript.GUIcallback(lstLuaFiles);

    /*
    glo_truss = makeTruss(truss);

    glo_capsula = glGenLists(1);
    glNewList( glo_capsula, GL_COMPILE );
    //Draw3D::drawCylinderStrip  ( 16, 10, 10, {0.0,0.0,-16.0,}, {0.0,0.0,16} );
    Draw3D::drawCapsula( (Vec3f){0.0,0.0,-1.0}, (Vec3f){0.0,0.0,1.0}, 2.0, 1.0, 0.7, 0.7, 0.2, 32, false );
    glEndList();
    //delete [] ups;
    */


    zoom = 1000.0;
}

//void SpaceCraftEditGUI::camera(){
//    camera_FreeLook( camPos );
//}

void SpaceCraftEditGUI::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    //glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//glDisable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_LIGHTING);
	glColor3f(1.0,0.5,1.0);
	if(glo_capsula) glCallList(glo_capsula);

	/*
    // to make nice antialiased lines without supersampling buffer
    // see  https://www.opengl.org/discussion_boards/showthread.php/176559-GL_LINE_SMOOTH-produces-bold-lines-%28big-width%29
    // https://www.opengl.org/discussion_boards/showthread.php/170072-GL_POLYGON_SMOOTH-is-it-that-bad
    if(bSmoothLines){
        glEnable (GL_BLEND);
        //glColor4f(0.0,0.0,0.0,0.1);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable (GL_LINE_SMOOTH);
        glHint (GL_LINE_SMOOTH_HINT, GL_FASTEST);
        glLineWidth(0.5);
        //glLineWidth(1.5);
	}
	//glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
    if(bWireframe){
        glEnable(GL_POLYGON_SMOOTH); // THIS WILL MAKE WIREFRAME AS SIDE EFFECT
        //glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);
        //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    }
    */

    if(glo_truss) glCallList(glo_truss);

    if(glo_ship) glCallList(glo_ship);

    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );

    //printf( "%i\n", EDIT_MODE::vertex );
    if(picked>=0){
        switch(edit_mode){
            case EDIT_MODE::vertex:
                glColor3f(1.0,1.0,1.0); Draw3D::drawPointCross( truss.points[picked], 0.3 );
                if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ){ Draw3D::drawLine(truss.points[picked], mouse_ray0); }
                break;
            case EDIT_MODE::edge  : glColor3f(1.0,1.0,1.0); auto ed = truss.edges[picked]; Draw3D::drawLine( truss.points[ed.a], truss.points[ed.b] ); break;
        }

    }

    mouse_ray0 = camMat.a*mouse_begin_x + camMat.b*mouse_begin_y;
    //glColor3f(0.0f,0.0f,0.0f); drawTruss( truss.edges.size(), &truss.edges[0], &truss.points[0] );
    //glColor3f(1.0f,1.0f,1.0f); Draw3D::drawPoints( truss.points.size(), &truss.points[0], 0.1 );

};

void SpaceCraftEditGUI::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

//void SpaceCraftEditGUI::keyStateHandling( const Uint8 *keys ){ };
/*
void SpaceCraftEditGUI::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    sprintf( curCaption, "%f %f\0", mouse_t, mouse_val );
    int ipoint_ = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
    if( (splines.ts[ipoint_+1]-mouse_t)<(mouse_t-splines.ts[ipoint_]) ) ipoint_++;
    if(ipoint_!=ipoint){
        ipoint=ipoint_;
        char buff[100];
        Vec3d r,v;
        r.set( splines.CPs[0][ipoint],splines.CPs[1][ipoint],splines.CPs[2][ipoint] );
        v.set( splines.getPointDeriv(ipoint,0), splines.getPointDeriv(ipoint,1), splines.getPointDeriv(ipoint,2) );
        sprintf(buff, "%i %f r(%3.3f,%3.3f,%3.3f) v(%3.3f,%3.3f,%3.3f)", ipoint, splines.ts[ipoint], r.x,r.y,r.z,   v.x,v.y,v.z );
        txtStatic.inputText = buff;
    }
    //printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
}
*/



void SpaceCraftEditGUI::keyStateHandling( const Uint8 *keys ){

    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ camPos.add_mul( camMat.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ camPos.add_mul( camMat.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ camPos.add_mul( camMat.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ camPos.add_mul( camMat.a, +0.05*zoom );  }

};

void SpaceCraftEditGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );

    gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                case SDLK_l:
                    //reloadShip( );
                    onSelectLuaShipScript.GUIcallback(lstLuaFiles);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    switch(edit_mode){
                        case EDIT_MODE::vertex: picked = truss.pickVertex( mouse_ray0, camMat.c, 0.5  ); printf("picked %i\n", picked); break;
                        case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    }; break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    switch(edit_mode){
                        case EDIT_MODE::vertex: int ip2 = truss.pickVertex( mouse_ray0, camMat.c, 0.5  ); if((picked>=0)&(ip2!=picked)); truss.edges.push_back((TrussEdge){picked,ip2,0}); break;
                        //case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    }; break;
                case SDL_BUTTON_RIGHT:break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

SpaceCraftEditGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);


	// https://www.opengl.org/discussion_boards/showthread.php/163904-MultiSampling-in-SDL
	//https://wiki.libsdl.org/SDL_GLattr
	//SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
    glEnable(GL_MULTISAMPLE);


	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new SpaceCraftEditGUI( junk , dm.w-150, dm.h-100 );
	//thisApp = new SpaceCraftEditGUI( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















