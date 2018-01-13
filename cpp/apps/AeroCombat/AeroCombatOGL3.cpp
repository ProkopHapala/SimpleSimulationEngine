

//#define SPEED_TEST


#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
#include <GL/glew.h>

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

#include "Terrain25D.h"
#include "Shooter.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "AeroCraftControl.h"
#include "AeroCraftWarrior.h"
//#include "AeroTest.h"

#include "DynamicControl.h"

//#include "AeroCraftControler.h"

//#include "FieldPatch.h"
#include "Solids.h"
//#include "AeroCraftWorld.h"
//#include "AeroCraftGUI.h"

/*
#include "SDL_utils.h"
//#include "Draw.h"
//#include "Draw3D.h"
//#include "AeroDraw.h"


#include "GUI.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "GUI.h"
#include "Plot2D.h"

#include "AeroCraftDesign.h"
*/


#include "Shader.h"
#include "GLObject.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "Shader.h"


// ====================================
//      AeroCraftGUI
// ====================================

class AeroCraftGUI : public AppSDL2OGL3, public SceneOGL3 { public:
	//AeroCraftWorld * world;
	//Shooter * world;

	Shooter            * world  = NULL;
    AeroCraftControler * pilot  = NULL;
    //AeroTester         * tester = NULL;

    DynamicControl rollControl;
    double roll;

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    /*
    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;
    GUIAbstractPanel*  focused = NULL;
    */

    bool mouseSteer   = false;
    bool autoRetractAirelon    = true;
    bool autoRetractRudder     = true;
    bool autoRetractElevator   = false;

    /*
    float  ElevatorRate = 0.01;
	float  AirelonRate  = 0.002;
	float  RudderRate   = 0.01;
    double AirelonRetractRate  = AirelonRate;
    double RudderRetractRate   = RudderRate;
    double ElevatorRetractRate = ElevatorRate;
    */

    bool useAutoPilot = false;

    float camDist = 100.0;

	int perFrame = 10;
	//double dt = 0.001;

	AeroCraftWarrior * myCraft;
	//AeroCraft * myCraft;
	AeroCraft * myCraft_bak;

	//FieldPatch fieldPatch;
	//int buildings_shape = -1;
	//int terrain_shape   = -1;

    bool staticTest = true;
    // - put to Spline manager ? ... make indepemented AeroCraft Test ?

	//int fontTex_DEBUG;

	AeroSurfaceDebugRecord leftWingRec,rightWingRec;


    Shader *sh1;
    GLMesh *glmesh,*gledges;
	/*
	Plot2D mainWingLD;
	Plot2D mainWingPolar;
	int polarPlotKind=1;

	QuePlot2D historyPlot;
	QuePlot2D wingsTrj;

	static const int nBaloons=100;
	Vec3d baloons[nBaloons];
	int gloBaloon;
	*/

	//static const int nProjectiles=100;
	//PointBody projectiles[nProjectiles];

	// ==== function declarations

	//void reallocateTrj(int n);
	//void resetSteer( );
	//void steerToDir( const Vec3d& dir );

	//void renderSkyBox( float x0, float y0, float z0 );
	//void drawStaticTest2D();

	//virtual void camera     ();
	//virtual void cameraHUD();

	//virtual void drawHUD();

	virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void draw( Camera& cam );

	AeroCraftGUI(int W, int H);

};

void AeroCraftGUI::update(){
    AppSDL2OGL3::update();

    Mat3d rot; rot.setT(myCraft->rotMat);
	world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
}

void AeroCraftGUI::draw( Camera& cam ){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //glEnable(GL_DEPTH_TEST);

    sh1->use();

    Mat3f mrot; mrot.setOne();
    sh1->set_modelMat( (GLfloat*)&mrot );
    sh1->set_modelPos( (const GLfloat[]){0.0f,0.0f,0.0f} );
    setCamera( *sh1, cam );

    int npoints = 1000;
    int nx = (int)(sqrt(npoints));
    int ny = npoints/nx;
    float step = 2.5;
    float span = 100.0;
    srand(15454);
    GLuint ucolor = sh1->getUloc("baseColor");
    glmesh->preDraw ();
    int nviewed = 0;
    for( int i=0; i<npoints; i++ ){
        //sh1->set_modelPos( (GLfloat*)(points+i) );
        //Vec3f pos = (Vec3f){ randf(-span,span),randf(-span,span),randf(-span,span) };
        //Vec3f pos = (Vec3f){ randf(-span,span),randf(-span,span),randf(-0,0) };
        Vec3f pos = (Vec3f){ step*(i/nx),step*(i%nx),0.0 };
        sh1->set_modelPos( (GLfloat*)&pos );
        //glUniform4f( ucolor, randf(0,1), randf(0,1), randf(0,1), 1.0 );
        //if( cam.pointInFrustrum(pos) )
        if( cam.sphereInFrustrum(pos,0.9) )
            { glUniform4f( ucolor, 1.0, 0.0, 0.0, 1.0 ); nviewed++; }
        else{ glUniform4f( ucolor, 0.0, 0.0, 1.0, 1.0 ); }
        printf( "nviewed %i \n", nviewed );
        //glUniform4fv( sh1->getUloc("baseColor"), 1,  (const float[]){1.0, 0.0, 0.0, 1.0} );
        //glmesh->draw();
        glmesh->drawRaw();
    }
};

AeroCraftGUI::AeroCraftGUI(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){
    /*
    printf("DEBUG 1 \n");

    initSDL( dm.w-150, dm.h-100);
    printf("DEBUG 2 \n");
    */

    for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );
    printf("DEBUG 3 \n");
    /*

    printf( " === GUI \n" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA.bmp" );
    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_Alpha.bmp" );
    panel.init( 5,5,105,35,  fontTex );
    panel.caption   = "rotation [Rad]"; panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;
    mpanel.initMulti( 120,5,200,120, fontTex , 4 );
    mpanel.caption="MultiPanel_1";
    txt.inputText = "insert number using =+-*";
    SDL_StartTextInput ();
    */

    printf( " === world  \n" );

	world = new Shooter();
    world->perFrame = 1;
    world->dt       = 0.005d;


    printf( " === aerocraft \n" );

    //char* fname = "data/AeroCraft1.ini";
    //char* fname = "data/AeroCraftStright1.ini";
    char* fname = "data/AeroCraftMainWing1.ini";
	myCraft_bak = new AeroCraft();          myCraft_bak->fromFile(fname);
	myCraft     = new AeroCraftWarrior();   myCraft    ->fromFile(fname);
	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");
    myCraft->pos.y=200.0;
    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );
    world->registrWarrior(myCraft);

    myCraft->leftAirelon->dbgRec  = &leftWingRec;
    myCraft->rightAirelon->dbgRec = &rightWingRec;

    printf( " === autoPilot1 \n" );

    pilot  = new AeroCraftControler();
    pilot->attach( myCraft );
    //pilot->craft0=myCraft_bak;

    pilot->rudder  .setSymetricRange(0.2);
    pilot->elevator.setSymetricRange(0.2);

    printf( " === tester \n" );

    /*
    tester      = new AeroTester();
    tester->autoPilot = pilot;
    tester->craft     = myCraft;
    tester->gravityG  = world->gravity;
    //tester->reallocateTrj(int n);
    */

    sh1=new Shader();
    sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    sh1->getDefaultUniformLocation();

    glmesh = new GLMesh();
    glmesh->init_wireframe( Solids::Octahedron );

    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0; cam.zmax = 1000.0; cam.zoom = 20.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;

};

void AeroCraftGUI::eventHandling( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            //case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            //case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;

            //case SDLK_LSHIFT   :
            //    //SDL_WarpMouseInWindow( window, WIDTH/2, HEIGHT/2);
            //    SDL_WarpMouseInWindow( window, WIDTH/2, HEIGHT*(1-pilot->elevator.getRelativeVal()) );
            //    mouseHandling();

            //case SDLK_u : useAutoPilot = !useAutoPilot;    break;
            //case SDLK_p : first_person = !first_person; break;
            //case SDLK_m : mouseSteer   = !mouseSteer;   break;
            //case SDLK_c : pilot->resetSteer( );
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    mouseSteer = true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    mouseSteer = false;
                    pilot->resetSteer();
                    break;
            }
            break;
    }


};

AeroCraftGUI * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new AeroCraftGUI( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}

/*
int main(int argc, char *argv[]){

    if (SDL_Init(SDL_INIT_VIDEO) < 0) die( "Unable to initialize SDL" );
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);
    window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if ( !window ) die("Unable to create window");
    context = SDL_GL_CreateContext( window );
    //SDL_GL_SetSwapInterval(1); // VSync On
    SDL_GL_SetSwapInterval(VSync);

    glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		quit();
		//return -1;
	}


	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	thisApp = new AeroCraftGUI( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->loop( 1000000 );
	return 0;
}
*/


