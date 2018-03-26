

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
#include "DrawOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "Shader.h"


#include "SimplexRuler.h"
#include "Ruler2DFast.h"
#include "TerrainHydraulics.h"

#include "TerrainOGL3.h"


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

    // Terrain - maybe move to Shooter
    SimplexRuler       ruler;
    //Ruler2DFast        square_ruler;
    TerrainHydraulics  hydraulics;
    double * ground    = NULL;

    TerrainOGL3 terrain1;


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


    Shader *sh1,*shDebug,*shTx;
    GLMesh *glmesh,*gledges,*msh_normals, *glDebug, *glTxDebug;
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

    shDebug=new Shader();
    //shDebug->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    shDebug->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shDebug->getDefaultUniformLocation();

    sh1=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    sh1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    sh1->getDefaultUniformLocation();

    shTx=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    shTx->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTx->getDefaultUniformLocation();

    sh1->use();
    glUniform3f(sh1->getUloc("lightColor"   ), 0.5,0.45,0.4 );
    glUniform3f(sh1->getUloc("diffuseColor" ), 1.0,1.0,1.0  );
    glUniform3f(sh1->getUloc("ambientColor" ), 0.2,0.25,0.3 );
    glUniform3f(sh1->getUloc("specularColor"), 2.0,2.0,2.0  );
    glUniform3f(sh1->getUloc("lightPos"     ), 10.0,-10.0,-10.0 );
    //glUniform3f(sh1->getUloc("lightPos"     ), 10.0,10.0,10.0 );

    GLMeshBuilder mshDebug;
    mshDebug.addLine      ( {0.0,0.0,0.0}, {10.0,10.0,10.0}, {1.0,0.0,0.0} );
    mshDebug.addPointCross( {0.0,0.0,0.0}, 1.0, {0.0,0.0,1.0} );
    glDebug = mshDebug.makeLineMesh();

    GLMeshBuilder mshbuild;

    //UVFunc2smooth( {10,10}, {0.0,-M_PI*0.0}, {M_PI*1.0,M_PI*0.5}, uvfunc , mshbuild );
    //UVFunc2wire( {10,10}, {0.0,-M_PI*0.5}, {M_PI*2.0,M_PI*0.5}, uvfunc , mshbuild );
    //UVFunc2wire( {10,10}, {0.0,-M_PI*0.5}, {M_PI*1.8,M_PI*0.4}, uvfunc , mshbuild );
    //UVFunc2smooth( {10,10}, {0.0,-M_PI*0.5}, {M_PI*1.8,M_PI*0.4}, uvfunc, mshbuild );
    //UVFunc2smooth( {20,20}, {0.0,-M_PI*0.499}, {M_PI*2.0,M_PI*0.499}, uvfunc, mshbuild );
    //UVFunc2smooth( {32,16}, {0.0,0.0}, {1.0,M_PI*2.00}, uvfunc, mshbuild );
    //UVFunc2smooth( {32,16}, {0.01,0.0}, {0.99,M_PI}, uvfunc, mshbuild );

    //Cone2Mesh( {20,20}, {0.01,0.0}, {0.99,M_PI*2.0},     1.0,0.2,3.0, false, mshbuild );
    //Sphere2Mesh  ( {20,20}, {-M_PI*0.4,0.0}, {M_PI*0.4,M_PI},     1.0        , false, mshbuild );
    //Torus2Mesh   ( {20,20}, {0.0,0.0},       {M_PI*2.0,M_PI*2.0}, 0.5,1.5,     false, mshbuild );
    Teardrop2Mesh( {20,16}, {0.01,0.0}, {0.99,M_PI*2.0},   0.8,0.2,6.0, false, mshbuild );    mshbuild.moveSub ( 0, {0.0,0.0,1.0} );
    //Sphere2Mesh  ( {20,20}, {-M_PI*0.49,0.0}, {M_PI*0.49,M_PI},     1.0 , false, mshbuild ); mshbuild.scaleSub( 1, {0.5,0.5,1.0} );
    Teardrop2Mesh( {10,8}, {0.01,0.0}, {0.99,M_PI}, 0.5,0.1,1.0, false, mshbuild ); mshbuild.moveSub( 1, {0.0,0.45,-1.5} );
    float naca1[4]={2.0,0.15,0.0,0.0};
    float naca2[4]={1.0,0.15,0.0,0.0};
    NACASegment2Mesh( {20,2}, {-1.0,0.0}, {1.0,1.0}, naca1,naca2, 5.0, false, mshbuild );
    //NACASegment2Mesh( {20,2}, {-1.0,0.0}, {1.0,1.0}, naca1,naca2, 5.0, false, mshbuild );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 3, {-1.0,1.0,1.0} );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 4, { 0.4,0.5,0.5} ); mshbuild.moveSub( 4, {0.0,0.0,-4.0} );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 5, {-0.4,0.5,0.5} ); mshbuild.moveSub( 5, {0.0,0.0,-4.0} );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 6, {-0.4,0.5,0.5} ); mshbuild.moveSub( 6, {0.0,0.0,-4.0} );  mshbuild.rotateSub( 6,{0.0,0.0,0.0},{0.0,0.0,1.0}, -M_PI*0.5 );

    //mshbuild.moveSub( 0, {1.0,2.0,3.0} );
    //mshbuild.rotateSub( 0, {0.0,0.0,0.0}, {1.0,0.0,0.0}, M_PI*0.5 );
    //mshbuild.scaleSub( 0, {1.0,0.5,0.25} );

    glmesh = mshbuild.makeGLmesh();
    //glmesh = mshbuild.normals2GLmesh(0.1);
    //glmesh->draw_mode = GL_LINES;
    msh_normals = mshbuild.normals2GLmesh(0.1);
    //glmesh = new GLMesh();
    //glmesh->init_wireframe( Solids::Octahedron );

    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0; cam.zmax = 1000.0; cam.zoom = 5.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;
    screen->camLookAt = new Vec3f(); (*screen->camLookAt)={0.0,0.0,0.0};


    double maxHeight = 1;

    int imgH = 1024;
    int imgW = 1024;

    ruler.setSize(imgW,imgH);
    ruler.setStep(50);
    //map_center = (Vec2d){ruler.na*0.75*ruler.step,ruler.nb*0.5*ruler.step};
    ground = new double[ruler.ntot];
    hydraulics.setSize(ruler.na,ruler.nb);
    hydraulics.ground = ground;
    //world.hydraulics.allocate( 512, 512 );
    hydraulics.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );

    for( int j=0; j<500; j++ ){
        int isz = 25;
        int ix0 = rand()%(hydraulics.nx-isz);
        int iy0 = rand()%(hydraulics.ny-isz);
        hydraulics.errodeDroples( 200, 100, 0.02, 0.15, 0.5, ix0, iy0, ix0+isz, iy0+isz );
    }
    for(int i=0; i<ruler.ntot; i++){ ground[i] *= maxHeight; };

    printf(" ruler.ntot %i \n", ruler.ntot );
    float* ground_f = new float[ruler.ntot];
    for(int i=0; i<ruler.ntot; i++){ ground_f[i] = (float)ground[i]; };


    /*
    int imgH = 1024;
    int imgW = 1024;
    float * height_map = new float[imgH*imgW];
    Vec2f step = (Vec2f){1.0/imgW,1.0/imgH};
    for( int iy=0; iy<imgH; iy++ ){
        for( int ix=0; ix<imgW; ix++ ){
            float x = ix*step.x-0.5;
            float y = iy*step.y-0.5;

            float r2 = x*x+y*y;
            float f = (1-4*r2);
            if(f<0){ f=0; }else{ f*=f; }
            //height_map[ iy*imgW + ix ] = (1/(1 + ))*sin(x*5)*sin(y*3)*0.5 + 0.5;

            f*= ( 1 + sin(x*5) * sin(y*8) +  cos(x*15) * cos(y*18)  )*0.3;

            height_map[ iy*imgW + ix ] = f;
        }
    }
    */
    //newTexture2D( txHeight, imgW, imgH, height_map, GL_RED, GL_FLOAT );
    //terrain1.init( {50,1000}, 1000.0,  {imgW,imgH},  height_map   );

    //terrain1.init( {200,1000}, 1000.0,  {imgW,imgH},  height_map   );
    terrain1.init( {200,1000}, 1000.0,  {imgW,imgH},  ground_f   );
    //glUniform3f( terrain1.ulocs.mapScale, 0.0, )
    terrain1.mapScale.z = 200.0;
    terrain1.mapScale.y = 0.0001;
    terrain1.mapScale.x = 0.0001;

    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);


    glTxDebug = new GLMesh();
    //DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts[]
    //glTxDebug->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs);
    glTxDebug->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs_2x2);

    //delete [] height_map;
    delete [] ground_f;


};

void AeroCraftGUI::update(){
    AppSDL2OGL3::update();

    Mat3d rot; rot.setT(myCraft->rotMat);
	//world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
}

void AeroCraftGUI::draw( Camera& cam ){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    sh1->use();
    Mat3f mrot; mrot.setOne();
    //sh1->set_modelMat( (GLfloat*)&mrot );
    //sh1->set_modelPos( (const GLfloat[]){0.0f,0.0f,0.0f} );

    //convert( myCraft->pos, cam.pos );
    cam.lookAt( myCraft->pos, 20.0 );
    setCamera( *sh1, cam );
    //setCameraOrtho(*sh1,cam);

    //sh1->setModelPose( myCraft->pos, myCraft->rotMat );
    sh1->setModelPoseT( myCraft->pos, myCraft->rotMat );

    //printf("orig "); myCraft->rotMat.printOrtho();

    //sh1->set_modelMat( myCraft->rotMat );
    //Mat3d r; r.setOne(); sh1->set_modelMat( r );
    //Mat3f r; r.setOne(); sh1->set_modelMat( (float*)&r );

    //sh1->set_modelMat( r );
    //sh1->setPose( myCraft->pos );

    GLuint ucolor = sh1->getUloc("baseColor");
    glUniform4f( ucolor, 0.0, 0.0, 0.0, 1.0 );  glmesh->draw();
    glUniform4f( ucolor, 1.0, 0.0, 0.0, 1.0 );  msh_normals->draw();

    shDebug->use();
    setCamera( *shDebug, cam );
    shDebug->setModelPoseT( myCraft->pos, myCraft->rotMat );
    glDebug->draw();

    //glmesh->draw(GL_TRIANGLES);

    //printf( "%f %f %f \n", cam.pos.x, cam.pos.y, cam.pos.z );

    // ==== terrain
    terrain1.pos.x = cam.pos.x;
    terrain1.pos.z = cam.pos.z;
    terrain1.pos.y = -100.0;
    terrain1.setViewRange( { cam.rot.c.x, cam.rot.c.z}, 0.3 );
    terrain1.sh.use();
    //terrain1.sh.set_camPos( (float*)&cam.pos );
    //terrain1.sh.set_camMat( (float*)&cam.rot );
    setCamera( terrain1.sh, cam );
    terrain1.draw();


    shTx->use();
    setCamera(*shTx, cam);
    shTx->setModelPoseT( myCraft->pos, myCraft->rotMat );

    /*
    float zoom = 100;
    Mat4f camMat; camMat.setOrthographic( zoom, zoom*cam.aspect, -1000, 1000 );
    shTx->set_camPos( (GLfloat*)&cam.pos );
    shTx->set_camMat( (GLfloat*)&camMat  );
    shTx->setModelPoseT( {0.0,0.0,0.0}, {1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0} );
    */


    glTxDebug->draw();



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

