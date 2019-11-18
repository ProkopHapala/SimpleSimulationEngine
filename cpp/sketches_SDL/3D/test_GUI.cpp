
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

#include "Table.h"

#include "GUI.h"
#include "Plot2D.h"

// ======================  TestApp

// References
// OpenGL GUI examples
//  * https://github.com/wjakob/nanogui
//  * https://github.com/ocornut/imgui

struct TestStruct{
    int    inum = 115;
    double dnum = 11.1154546;
    double fvoid= 0.0;
    float  fnum = 11.115;
    Vec3d  dvec = (Vec3d){ 1.1545, 2.166, 3.1545};
};

void command_example( double value ){
    printf( "I'm writting value : %6.6f \n", value );
}

class TestAppGUI : public AppSDL2OGL_3D {
	public:
    int      fontTex;
    int      fontTex3D;

    GUI gui;
    //GUIPanel   panel;
    //MultiPanel mpanel;
    ScisorBox clipBox;
    //DropDownList list1;
    //GUIAbstractPanel*  focused = NULL;

    GUITextInput txt;
    Plot2D plot1;

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling ( const SDL_Event& event  );

	TestAppGUI( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppGUI::TestAppGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    GUI_fontTex = fontTex;
    fontTex3D = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ====== Table
    Table* tab1 = new Table();
    tab1->n      =  120;
    TestStruct* tab_data = new TestStruct[tab1->n];
    for(int i=0; i<tab1->n; i++){ tab_data[i].inum = i; tab_data[i].fnum = i*0.1; tab_data[i].dnum = sq(i*0.1); }
    //char* tab_data_ = (char*)tab_data;
    //char* pinum     = (char*)&(tab_data->inum);
    //int ipinum      = pinum - tab_data_;
    //tab1->data       = (char*)tab_data;
    tab1->bind(tab_data, sizeof(*tab_data) );
    /*
    tab1->columns.push_back( Atribute( ((char*)&(tab_data->inum))-((char*)tab_data), 1, DataType::Int   ) );
    tab1->columns.push_back( Atribute( ((char*)&(tab_data->fnum))-((char*)tab_data), 1, DataType::Float ) );
    tab1->columns.push_back( Atribute( ((char*)&(tab_data->dnum))-((char*)tab_data), 1, DataType::Double) );
    tab1->columns.push_back( Atribute( ((char*)&(tab_data->dvec))-((char*)tab_data), 3, DataType::Double) );
    */

    tab1->addColum( &(tab_data->inum), 1, DataType::Int    );
    tab1->addColum( &(tab_data->fnum), 1, DataType::Float  );
    tab1->addColum( &(tab_data->dnum), 1, DataType::Double );
    tab1->addColum( &(tab_data->dvec), 3, DataType::Double );

    gui.addPanel( new TableView( tab1, "tab1", 150.0, 250.0,  0, 0, 5, 3 ) );

    return;

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

	Draw3D::drawAxis ( 3.0f );

	glColor4f(1.0f,1.0f,1.0f,0.9f);
	txt.view3D( {1.0,1.0,1.0}, fontTex3D, 0.3 );

	clipBox.apply();
	Draw2D::drawText("jkfbksfdfjshbdfs \n sdjfhsdjkfksdh \n djfshdkjfhsdkhf", (Vec2d){-6.0,5.0},(Vec2d){12.0,12.0}, fontTex3D, 1.0 );
    glDisable(GL_SCISSOR_TEST);

};

void TestAppGUI::drawHUD(){

    glColor3f(1.0,1.0,1.0);
    txt.viewHUD( {100,220}, fontTex );

	gui.draw();

	glTranslatef( 10.0,300.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, fontSizeDef, {10,5} );

	glTranslatef( 200.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	plot1.view();

}

void TestAppGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
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
















