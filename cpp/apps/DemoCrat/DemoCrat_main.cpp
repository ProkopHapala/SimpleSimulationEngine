
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <dlfcn.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"
#include "SDL_utils.h"
#include "runtime.h"
#include "Demo.h"

#include "GUI.h"
#include "Plot2D.h"
#include "GUI.h"
#include "IO_utils.h"

// see
//   http://www.yolinux.com/TUTORIALS/LibraryArchives-StaticAndDynamic.html
// check symbols
//   readelf -s /usr/lib64/libjpeg.so
//   readelf -s libDemo1.so | grep draw
//   readelf -s libDemo1.so | grep setup
// must be compiled with
//   g++ -rdynamic $(FLAGS) -o program main.cpp -ldl -lm -lSDL2 -lGL -lGLU
//   see:  /home/prokop/Dropbox/MyDevSW/Ccko/RuntimePluginCPP
//   see: https://stackoverflow.com/questions/11783932/how-to-add-linker-or-compile-flag-in-cmake-file

//int  default_font_texture;
char str[0x10000];

//const char * fflags = "-I.. -I/home/prokop/git/SimpleSimulationEngine/cpp/common/math ";
//const char * fflags = "-Wall -g -Og  -I.. -I../../../../common/math -I../../../../common/utils ";
const char * fflags = "-w -I.. -I../../../../common/math -I../../../../common/utils ";

class DemoCratApp : public AppSDL2OGL { public:

    void        *lib_handle = 0;
    Pprocedure   pplSetup = 0;
    Pprocedure   pplDraw = 0;
    PmouseFunc   pplOnMouse = 0;
    PDemoFactory pDemoFactory = 0;
    Demo*        demo = 0;

    class OnSelectDemo : public GUIEventCallback{ public:
        DemoCratApp* app;
        //OnSelectDemo(DemoCratApp* app_){app=app_;}
        virtual int GUIcallback(GUIAbstractPanel* caller) override {
            //((DropDownList*)caller)->selectedToStr(str+sprintf(str,"data/"));
            ((DropDownList*)caller)->selectedToStr(str);
            app->loadDemo(str);
        }
    };

    GUI gui;
    DropDownList* lstDemos;
    OnSelectDemo  onSelectDemo;

    //int bRecompile = false;
    int bRecompile = true; // Force recompile

    int loadDemo( char * fname );

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling( const SDL_Event& event  );
	void debug_buffinsert( );
	//void pickParticle( Particle2D*& picked );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

	double polyLinIntercept( int n, double* xs, double* ys, double y0, double a );
	void   drawVisibilityIsolines( Vec2d ray0, int ndhs, int nDirs, double phiMin, double phiMax, double dhmin, double dhmax, double tmax );

	DemoCratApp( int& id, int WIDTH_, int HEIGHT_ );

};

int DemoCratApp::loadDemo( char * fname ){
    printf( "DemoCratApp::loadDemo '%s'", fname );
    int slen = strlen(fname); if(slen>4) fname[slen-4] = '\0';
    printf( "DemoCratApp::loadDemo '%s'", fname );

    char *error;
    demo = 0;
    pDemoFactory=0;
    //pplDrawXY=0;
    pplSetup=0;
    pplDraw=0;

    if(lib_handle){ dlclose(lib_handle); lib_handle=0; };

    char fullname[1024];
    sprintf(fullname, "data/bin/lib%s.so", fname );

    if( bRecompile || !fileExist(fullname) ) recompileLib( "data/bin", "../", fname, fflags );

    printf( "=== readelf -s data/bin/libDemo1.so | grep draw\n" );
    system( "readelf -s data/bin/libDemo1.so | grep draw" );
    printf( "=== readelf -s data/bin/libDemo1.so | grep setup\n" );
    system( "readelf -s data/bin/libDemo1.so | grep setup" );

    printf("loading %s :\n", fullname );

    void* plib = dlopen( fullname, RTLD_LAZY);
    //if (!plib){ printf( "%s\n", dlerror()); exit(1); }
    if (plib){
        //if(lib_handle){
        //    printf( "dlclose %i \n", lib_handle );
        //    dlclose(lib_handle);
        //};
        lib_handle=plib;
    }else{ printf( "%s\n", dlerror()); return -1; }

    pplSetup = (Pprocedure)dlsym(lib_handle, "plSetup");
    if ((error = dlerror())){ printf("%s\n", error); pplSetup=0; }

    pplDraw = (Pprocedure)dlsym(lib_handle, "plDraw");
    if ((error = dlerror())){ printf( "%s\n", error); pplDraw=0; }

    pplOnMouse = (PmouseFunc)dlsym(lib_handle, "plOnMouse");
    if ((error = dlerror())){ printf( "%s\n", error); pplOnMouse=0; }

    pDemoFactory = (PDemoFactory)dlsym(lib_handle, "CreateDemo");
    if ((error = dlerror())){ printf( "%s\n", error); pDemoFactory=0; }

    if(pDemoFactory){
        demo = pDemoFactory();
        demo->setup();
    };

    if(pplSetup) pplSetup();
    return 0;
}


DemoCratApp::DemoCratApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    //default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //onSelectDemo = new OnSelectDemo(this);
    lstDemos  = new DropDownList( "lua files",20,HEIGHT_-100,200,5);
    gui.addPanel(lstDemos);
    listDirContaining( "data", ".cpp", lstDemos->labels );
    onSelectDemo.app = this; lstDemos->onSelect = &onSelectDemo;

    //printf( "default_font_texture :  %i \n", default_font_texture );
}

void DemoCratApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	if(demo) demo->onMouse(mouse_begin_x,mouse_begin_y,0);
	if(demo) demo->draw();
};

void DemoCratApp::drawHUD(){
    gui.draw();
};

void DemoCratApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    gui.onEvent(mouseX,HEIGHT-mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                case SDLK_l:
                    //reloadShip( );
                    onSelectDemo.GUIcallback(lstDemos);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    break;
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
    camStep = zoom*0.05;
}

// ===================== MAIN

DemoCratApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);

    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    //printf( "Display size W: %i H:%i \n", dm.w, dm.h );

	int junk;
	thisApp = new DemoCratApp( junk , dm.w-150, dm.h-100 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















