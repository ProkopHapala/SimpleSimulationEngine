
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec2.h"
#include "Solids.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "GUI.h"
#include "Plot2D.h"

// ======================  TestApp

// References
// OpenGL GUI examples
//  * https://github.com/wjakob/nanogui
//  * https://github.com/ocornut/imgui

void command_example( double value ){
    printf( "I'm writting value : %6.6f \n", value );
}

class TestAppGUI : public AppSDL2OGL_3D {
	public:
    int      fontTex;
    int      fontTex3D;


    GUIPanel   panel;
    MultiPanel mpanel;
    ScisorBox clipBox;
    DropDownList list1;
    GUIAbstractPanel*  focused = NULL;

    GUI gui;


    GUITextInput txt;
    Plot2D plot1;

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling ( const SDL_Event& event  );

	TestAppGUI( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppGUI::TestAppGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    GUI_fontTex = fontTex;
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_Alpha.bmp" );
    fontTex3D = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    /*
    list1.initList("DropList1",100.0,300.0,200.0,5);
    list1.labels.push_back("Item_1");
    list1.labels.push_back("Item_2");
    list1.labels.push_back("Item_3");
    list1.labels.push_back("Item_4");
    list1.labels.push_back("Item_5");
    list1.labels.push_back("Item_7");
    list1.labels.push_back("Item_8");
    //list1.caption = "DropList1";
    gui.addPanel( &list1 );
    */
    ((DropDownList*)gui.addPanel( new DropDownList("DropList1",100.0,300.0,200.0,5) ))
        ->addItem("Item_1")
        ->addItem("Item_2")
        ->addItem("Item_3")
        ->addItem("Item_4")
        ->addItem("Item_5")
        ->addItem("Item_6")
        ->addItem("Item_7")
        ->addItem("Item_8");



    clipBox.initScisor("Scissor_1",300,200,500,400);
    //clipBox.caption = "Scissor_1";
    gui.addPanel( &clipBox );
    //gui.addPanel( new ScisorBox("Scissor_1",300,200,500,400) );

    /*
    panel.initPanel( "rotation [Rad]", 5,5,105,35 );
    //panel.caption   = "rotation [Rad]";
    //panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;
    //sprintf(panel.val_text, "NA" );
    panel.command = &command_example; panel.isButton = true;
    gui.addPanel( &panel );
    */
    ((GUIPanel*)gui.addPanel( new GUIPanel( "rotation [Rad]", 5,5,105,35, true, true ) ))
        ->command = &command_example;

    /*
    //mpanel.initMulti( 120,5,200,120, fontTex , 4 );
    mpanel.initMulti( "MultiPanel_1", 120,5,200,fontSizeDef*4, 4 );
    //mpanel.caption="MultiPanel_1";
    mpanel.bgColor = 0x9090A0;
    gui.addPanel( &mpanel );
    */
    gui.addPanel( new MultiPanel( "MultiPanel_1", 120,5,200,fontSizeDef*4, 4 ) )
        ->bgColor = 0x9090A0;

    txt.inputText = "insert number using =+-*/";

    SDL_StartTextInput ();
    //panel.nChars = 6;

    // Plot2D

    plot1.init();
    plot1.fontTex = fontTex;
    plot1.clrGrid = 0xFF404040;
    plot1.clrBg   = 0xFF408080;

    DataLine2D * line1 = new DataLine2D(100);
    line1->linspan(-3*M_PI,2*M_PI);
    //Func1d myFunc = &sin;
    //line1->yfunc = myFunc;
    line1->yfunc = &sin;

    DataLine2D * line2 = new DataLine2D(400);
    line2->linspan(-M_PI,3*M_PI);
    for(int i=0; i<line2->n; i++){
        double x     = line2->xs[i];
        line2->ys[i] = sin(x*x);
    }
    line2->clr = 0xFF00FF00;

    plot1.lines.push_back( line1 );
    plot1.lines.push_back( line2 );
    plot1.render();

}

void TestAppGUI::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable ( GL_LIGHTING );

	//Draw3D::drawText( txt.inputText.c_str(), {1.0,1.0,1.0}, fontTex, 0.3, 0, 0 );

	Draw3D::drawAxis ( 3.0f );

	glColor4f(1.0f,1.0f,1.0f,0.9f);
	txt.view3D( {1.0,1.0,1.0}, fontTex3D, 0.3 );

	//Rect2d rect = (Rect2f){-5.0,-5,5,5};
	//Vec2f p1=(Vec2f){-5.0,-5.0},p2=(Vec2f){5.0,5.0};
	//glColor3f(1.0,0.0,1.0);
	//Draw2D::drawRectangle(p1,p2, false );

	//glEnable(GL_SCISSOR_BOX);
	//glEnable(GL_SCISSOR_TEST);
	//glScissor(300,  200,  200,  200);
	clipBox.apply();
	//Draw3D::drawSphere_oct(5,4.0, (Vec3d){0.0,0.0,0.0} );
	Draw2D::drawText("jkfbksfdfjshbdfs \n sdjfhsdjkfksdh \n djfshdkjfhsdkhf", (Vec2d){-6.0,5.0},(Vec2d){12.0,12.0}, fontTex3D, 1.0 );
    glDisable(GL_SCISSOR_TEST);
	//Draw3D::drawText( "AHOJ!\0", {1.0,1.0,1.0}, fontTex, 5.0, 0, 0 );

};

void TestAppGUI::drawHUD(){

    glColor3f(1.0,1.0,1.0);
    txt.viewHUD( {100,220}, fontTex );


    /*
    //Draw2D::drawRectangle( {300.0,200.0},{500.0,400.0}, false );
    //panel.tryRender(); clipBox.draw();
    clipBox.draw();
	//panel .tryRender();  panel.draw();
	//mpanel.tryRender(); mpanel.draw();
	list1.draw();
    panel.draw();
	mpanel.draw();
	if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);
	*/
	gui.draw();

	glTranslatef( 10.0,300.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, fontSizeDef, {10,5} );

	glTranslatef( 200.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	plot1.view();

	//glColor3f(1.0f,1.0f,1.0f);
	//panel.render();
	//Draw2D::drawRectangle ( 100, 100, 300, 200, false    );
	//Draw2D::drawLine      ( {0.0f,0.0f},{100.0f, 200.0f} );
}

void TestAppGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    /*
    switch( event.type ){
        case SDL_KEYDOWN:
            if(focused){ focused->onKeyDown( event ); }else{ txt.onKeyDown(  event ); }; break;
        case SDL_TEXTINPUT:
            if(focused){ focused->onText   ( event ); }else{ txt.onText   ( event );  }; break;
        case SDL_MOUSEBUTTONDOWN:
            //printf( "%i %i\n", mouseX, mouseY );
            GUIAbstractPanel* active = NULL; focused=NULL;
            active = mpanel .onMouse( mouseX, mouseY, event ); if(active) focused=active;
            active = panel  .onMouse( mouseX, mouseY, event ); if(active) focused=active;
            active = clipBox.onMouse( mouseX, mouseY, event ); if(active) focused=active;
            active = list1  .onMouse( mouseX, mouseY, event ); if(active) focused=active;
            break;
    };
    */
    gui.onEvent( mouseX, mouseY, event );
    //GUI_globalEventHandler( &event );
    AppSDL2OGL::eventHandling( event );
    //camStep = zoom*0.05;
}

// ===================== MAIN

TestAppGUI * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppGUI( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















