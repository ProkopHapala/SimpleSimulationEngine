

//#define SPEED_TEST

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Draw3D.h"
//#include "testUtils.h"

#include "Body.h"

//#include "SimplexRuler.h"
//#include "Ruler2DFast.h"
//#include "TerrainHydraulics.h"

//#include "Terrain25D.h"
//#include "Shooter.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "AeroCraftControl.h"
#include "AeroCraftWarrior.h"

//#include "AeroCraftControler.h"

//#include "FieldPatch.h"
#include "Solids.h"
//#include "AeroCraftWorld.h"
//#include "AeroCraftEditor.h"

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AeroDraw.h"
#include "AeroTest.h"

#include "GUI.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "GUI.h"
#include "Plot2D.h"

#include "AeroCraftDesign.h"

// ===============================
// ===== GLOBAL CONSTAMNTS
// ===============================

// ===============================
// ===== Free Functions
// ===============================

// ====================================
//      AeroCraftEditor
// ====================================

class AeroCraftEditor : public AppSDL2OGL_3D { public:

    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;
    GUIAbstractPanel*  focused = NULL;

	AeroCraftWarrior * myCraft;

    const Uint8 *scanKeys;
    Uint32 mouseButtons;
    Vec3d mouseRay0;

    bool staticTest = true;
    // - put to Spline manager ? ... make indepemented AeroCraft Test ?

	//int fontTex_DEBUG;

	AeroSurfaceDebugRecord leftWingRec,rightWingRec;
	Plot2D mainWingLD;
	Plot2D mainWingPolar;
	int polarPlotKind=1;

	// ==== function declarations

	//virtual void camera     ();
	//virtual void cameraHUD();
	virtual void draw   ();
	virtual void drawHUD();

    virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling   ( );

	AeroCraftEditor( int& id, int WIDTH_, int HEIGHT_ );

};

void AeroCraftEditor::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    //glShadeModel(GL_FLAT);
    glShadeModel(GL_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glColor3f(1.0,1.0,1.0);

    renderAeroCraft(*myCraft, false);

    mouseRay0 = camPos + camMat.a*mouse_begin_x + camMat.b*mouse_begin_y;
    //printf( "mouse_begin_x" );
    Draw3D::drawPointCross( mouseRay0, 1.0 );

};

void AeroCraftEditor::drawHUD(){
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);


	//panel .tryRender();  panel.draw();
	//mpanel.tryRender(); mpanel.draw();
	//if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);

    float sc = WIDTH*0.04;
    glPushMatrix();
    if(polarPlotKind==1){
        glTranslatef(200.0,HEIGHT-150.0,200.0);
        glScalef    (sc*4,sc,WIDTH);
        mainWingPolar.view();
        Draw2D::drawLine( {0.0,-1.0}, {0.0,1.0} );
        glColor3f(1.0,1.0,0.0); Draw2D::drawPointCross( {leftWingRec.CD, leftWingRec.CL}, 0.05 );
        glColor3f(0.0,1.0,1.0); Draw2D::drawPointCross( {rightWingRec.CD,rightWingRec.CL}, 0.05 );
    }else if( polarPlotKind==2 ){
        double phiL = atan2( leftWingRec.sa,  leftWingRec.ca );
        double phiR = atan2( rightWingRec.sa, rightWingRec.ca );
        glTranslatef(200.0,HEIGHT-150.0,200.0);
        glScalef    (sc,sc,WIDTH);
        mainWingLD.view();
        glColor3f(1.0,1.0,0.0); mainWingLD.drawVline(phiL);
        glColor3f(0.0,1.0,1.0); mainWingLD.drawVline(phiR);
    }
    glPopMatrix();

}

AeroCraftEditor:: AeroCraftEditor( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA.bmp" );
    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_Alpha.bmp" );
    panel.init( 5,5,105,35,  fontTex );
    panel.caption   = "rotation [Rad]"; panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;
    mpanel.initMulti( 120,5,200,120, fontTex , 4 );
    mpanel.caption="MultiPanel_1";
    txt.inputText = "insert number using =+-*/";
    SDL_StartTextInput ();

    printf( " === aerocraft \n" );

    //char* fname = "data/AeroCraft1.ini";
    char* fname = "data/AeroCraftStright1.ini";
	myCraft     = new AeroCraftWarrior();   myCraft    ->fromFile(fname);
	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");

    printf( " DEBUG 1 \n" );
    // Polar Plotting
    mainWingLD.init();
    mainWingLD.fontTex = fontTex;
    mainWingLD.clrGrid = 0xFF404040;
    //mainWingLD.clrBg   = 0xFF408080;
    printf( " DEBUG 2 \n" );
    int nsamp = 100;
    double phiRange = M_PI*0.5;
    DataLine2D * lLift = new DataLine2D(nsamp); mainWingLD.lines.push_back( lLift ); lLift->linspan(-phiRange,phiRange); lLift->clr = 0xFFff0000;
    DataLine2D * lDrag = new DataLine2D(nsamp); mainWingLD.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;
    printf( " DEBUG 3 \n" );
    mainWingPolar.init();
    mainWingPolar.fontTex = fontTex;
    mainWingPolar.clrGrid = 0xFF404040;
    DataLine2D * LDpolar = new DataLine2D(nsamp); mainWingPolar.lines.push_back( LDpolar );
    printf( " DEBUG 4 \n" );
    for(int i=0; i<nsamp; i++){
        double phi = lLift->xs[i];
        double ca = cos(phi);
        double sa = sin(phi);
        double CD,CL;
        myCraft->leftAirelon->polarModel(ca,sa, CD, CL);
        lLift->ys[i]   = CL;
        lDrag->ys[i]   = CD;
        LDpolar->xs[i] =CD;
        LDpolar->ys[i] =CL;

    }
    printf( " DEBUG 5 \n" );
    mainWingLD.update();
    mainWingLD.autoAxes(0.5,0.2);
    mainWingLD.render();
    printf( " DEBUG 6 \n" );
    mainWingPolar.update();                 printf( " DEBUG 6.1 \n" );
    mainWingPolar.autoAxes(0.1,0.5);        printf( " DEBUG 6.2 \n" );
    mainWingPolar.render();                 printf( " DEBUG 6.3 \n" );
    printf( " DEBUG 7 \n" );
    printf( " INITIALIZATION DONE \n" );
};

void AeroCraftEditor:: eventHandling   ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
    }
};

void AeroCraftEditor:: keyStateHandling( const Uint8 *keys ){
    scanKeys = keys;
    AppSDL2OGL_3D::keyStateHandling( keys );
};

AeroCraftEditor * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	thisApp = new AeroCraftEditor( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->loop( 1000000 );
	return 0;
}


