
//#define SPEED_TEST

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
#include <GL/glew.h>

#include "testUtils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "VecN.h"

//#include "Body.h"
//#include "DynamicControl.h"
//#include "FieldPatch.h"
#include "Solids.h"

#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "Shader.h"
#include "DrawOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"


#define _DEBUG_VIEW_ 1
#include "DebugView.h"  //do we need it ?
DEBUG_VIEW_DEFINE()

#include "SDL_utils.h"
#include "IO_utils.h"


#include "MusicUtils.h"
#include "MusicRendererOGL3.h"



using namespace Music;


// ====================================
//      MusicVisualizerGUI
// ====================================

class MusicVisualizerGUI : public AppSDL2OGL3, public SceneOGL3 { public:

    //DynamicControl rollControl;

    Spectrum waveform;

    float dtJul = 0.2;
    Vec2f CJul = Vec2fZero;

    float camDist = 100.0;
	int perFrame = 10;
	//double dt = 0.001;

	//int fontTex_DEBUG;

    ShaderStack layers;


    Shader *sh1=0,*shDebug=0,*shTx=0,*shJulia=0,*shTex=0;
    GLMesh *histMesh=0, *waveMesh=0, *glmesh=0,*gledges=0,*msh_normals=0, *glDebug=0, *glTxDebug=0;


	// ==== function declarations

	//void reallocateTrj(int n);

	//void renderSkyBox( float x0, float y0, float z0 );

	//virtual void camera     ();
	//virtual void cameraHUD();

	//virtual void drawHUD();

	//virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );

    virtual void draw( Camera& cam );

	MusicVisualizerGUI(int W, int H);

};

MusicVisualizerGUI::MusicVisualizerGUI(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){

 for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    DEBUG_VIEW_INIT()


    //frameBuff1.init( WIDTH, HEIGHT );
    //texTest = makeTestTextureRGBA( 256, 256);

    shTex=new Shader();
    shTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTex->getDefaultUniformLocation();



    //for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    shDebug=new Shader();
    shDebug->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    //shDebug->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shDebug->getDefaultUniformLocation();

    sh1=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    sh1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    sh1->getDefaultUniformLocation();

    /*
    shTx=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    shTx->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTx->getDefaultUniformLocation();
    */

    shJulia=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    //shJulia->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shJulia->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/Julia.glslf"   );
    shJulia->getDefaultUniformLocation();

    sh1->use();
    glUniform3f(sh1->getUloc("lightColor"   ), 0.5,0.45,0.4 );
    glUniform3f(sh1->getUloc("diffuseColor" ), 1.0,1.0,1.0  );
    glUniform3f(sh1->getUloc("ambientColor" ), 0.2,0.25,0.3 );
    glUniform3f(sh1->getUloc("specularColor"), 2.0,2.0,2.0  );
    glUniform3f(sh1->getUloc("lightPos"     ), 10.0,-10.0,-10.0 );
    //glUniform3f(sh1->getUloc("lightPos"     ), 10.0,10.0,10.0 );

    GLMeshBuilder mshDebug;
    mshDebug.addLine      ( (Vec3f){0.0,0.0,0.0}, {10.0,10.0,10.0}, {1.0,0.0,0.0} );
    mshDebug.addPointCross( {0.0,0.0,0.0}, 1.0, {0.0,0.0,1.0} );
    glDebug = mshDebug.makeLineMesh();

    // ========== Mesh Builder

    GLMeshBuilder mshbuild;

    Parabola2Mesh ( {20,10}, {-1.5,-0.5*M_PI}, {1.0,0},      2.0, 4.0, 0.0, false, mshbuild );
    Hyperbola2Mesh( {20,20}, {-1.5,0.0}, {1.0,M_PI},    1.5, 2.0, 4.0, 0.0, false, mshbuild );

    //mshbuild.moveSub( 0, {1.0,2.0,3.0} );
    //mshbuild.rotateSub( 0, {0.0,0.0,0.0}, {1.0,0.0,0.0}, M_PI*0.5 );
    //mshbuild.scaleSub( 0, {1.0,0.5,0.25} );

    glmesh      = mshbuild.makeGLmesh();
    msh_normals = mshbuild.normals2GLmesh(0.1);



    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0;
    cam.zmax = 1000.0;
    cam.zoom = 5.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;



    glTxDebug = new GLMesh();
    //DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts[]
    //glTxDebug->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs);
    glTxDebug->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs_2x2);

    layers.screenQuad = glTxDebug;
    layers.makeBuffers( 2, screen[0].WIDTH, screen[0].HEIGHT );
    layers.shaders.push_back( shTex   );
    layers.shaders.push_back( shJulia );

    //return 0;


    int result = 0;
    int flags = MIX_INIT_MP3;

    const char *music_file_name = "common_resources/02-Lazenca-SaveUs.mp3";
    //const char *music_file_name = "Yanni - Reflections Of Passion.mp3";

    Mix_OpenAudio(22050, AUDIO_S16SYS, 2, 640);
    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        printf("Failed to init SDL\n");
        exit(1);
    }
    if (flags != (result = Mix_Init(flags ))) {
        printf("Could not initialize mixer (result: %d).\n", result);
        printf("Mix_Init: %s\n", Mix_GetError());
        exit(1);
    }
    Mix_Music *music = Mix_LoadMUS( music_file_name );


    waveform.realloc( 1024, 16 );
    waveform.mixer_info();
    waveform.clearHist();

    Vec3f ps[3*waveform.nwave];
    //for(int i=0; i<waveform.nwave; i++){ ps[i].set( i*0.1, waveform.wave[i]*0.0001, 0 ); }
    for(int i=0; i<waveform.nwave; i++){ ps[i].set( 0.0, 0, 0 ); }
    waveMesh = polyLineMesh( waveform.nwave, (float*)ps );
    histMesh = polyLineMesh( waveform.nhist, (float*)ps );

    Mix_SetPostMix( postmix_Spectrum, (void*)&waveform );
    Mix_PlayMusic(music, 1);


    //camStep = 2.0;
};


void MusicVisualizerGUI::draw( Camera& cam ){

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    //waveform.spectrumHistSmearing();
    //waveform.update( 5.1 );
    waveform.need_refresh = false;

    //fflush( stdout );
    printf("\e[1;1H\e[2J"); // clear screen //https://stackoverflow.com/questions/2347770/how-do-you-clear-the-console-screen-in-c
    waveform.printSpectrum();
    //printf( "nwave %i nhist %i \n", waveform.nwave, waveform.nhist );


    shDebug->use();
    cam.lookAt( (Vec3f){0.0,0.0,0.0}, 20.0 );

    setCamera( *shDebug, cam );
    shDebug->setModelPoseT( (Vec3d){-6.0f,0.0,0.0}, Mat3dIdentity );


    //plotBuffStereo ( *waveMesh, *shDebug, waveform.nwave, waveform.wave, 0.05, 0.0001 );

    //FFT(  waveform.wave, waveform.nwave, 1 );
    Vec2d* wave = (Vec2d*)(waveform.wave);
    //for(int i=0; i<waveform.nwave+16; i++){ wave[i].fromAngle( i*0.04f); wave[i].add( Vec2d::newFromAngle( i*0.1f) ); wave[i].mul(1e+6); };
    //FFT(  waveform.wave+2, waveform.nwave/2, 1 );
    //plotBuffStereo ( *waveMesh, *shDebug, waveform.nwave, waveform.wave, 0.01, 0.0001 );
    //FFT(  waveform.wave+2, waveform.nwave/2, 1 );


    float dx = 0.01;
    waveform.update( 5.1 );
    plotBuffStereo ( *waveMesh, *shDebug, waveform.nwave, waveform.Fwave, dx, 0.000001 );


    //fflush( stdout );
    printf("\e[1;1H\e[2J"); // clear screen //https://stackoverflow.com/questions/2347770/how-do-you-clear-the-console-screen-in-c
    waveform.printSpectrum();
    shDebug->setUniformVec4f("baseColor", (Quat4f){0.0,1.0,0.0,1.0});
    plotBuff( *histMesh, waveform.nhist, waveform.hist, 0.5*dx*waveform.nwave/waveform.nhist, 0.000001 );


    //plotBuffStereo ( *waveMesh, *shDebug, 396, waveform.wave, 0.05, 0.00001 );
    //VecN::set( waveform.nwave*2, 0.0, waveform.wave );

    int narg;



    //  Following will be rendered to BUFFER[1]
    layers.bindOutput( 0 );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);


    shJulia->use();
    //uint locC = shJulia->getUloc("C");
    float scJulC = 0.0000005;
    //shJulia->setUniformVec2f( "Const", {waveform.hist[0]*scJulC,waveform.hist[waveform.nwave/2]*scJulC} );
    float Cx=2-waveform.hist[0               ]*scJulC*0.2;
    float Cy=  waveform.hist[waveform.nhist/2]*scJulC;

    CJul.mul(1-dtJul); CJul.add(Cx*dtJul,Cy*dtJul);
    printf("C %g %g | %g %g \n", CJul.y, CJul.y, Cx,Cy );
    //shJulia->setUniformVec2f( "Const", {Cx,Cy} );
    shJulia->setUniformVec2f( "Const", CJul );
    //shJulia->setUniformVec2f( "Const", {sin(frameCount*0.02),cos(frameCount*0.01)} );
    //shJulia->setUniformVec2f( "Const", {-0.3,0.6} );
    setCamera(*shJulia, cam);
    shJulia->setModelPoseT( (Vec3d){-4.,-4.,0.0}, Mat3dIdentity*8.0 );
    glTxDebug->draw();


    layers.unbindOutput();



    /*
    shTex->use();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture  (GL_TEXTURE_2D, layers.buffers[0]->texRGB   );
    setCamera(*shTex, cam );
    shTex->setModelPoseT( (Vec3d){-4.,-4.,0.0}, Mat3dIdentity*8.0 );
    glTxDebug->draw();
    */


    shTex->use();
    setCamera(*shTex, cam );
    shTex->setModelPoseT( (Vec3d){-4.,-4.,0.0}, Mat3dIdentity*8.0 );
    layers.render( 0, -1, 1, (const int[]){0} ); // renders using shader[0] to default_screen using 1 input buffer #0




};


/*
void MusicVisualizerGUI::update(){
    AppSDL2OGL3::update();

    //mouseButtons = SDL_GetMouseState(&mx,&my);
    bool RMB = mouseButtons&SDL_BUTTON(SDL_BUTTON_RIGHT);

    if      ( keys[ SDL_SCANCODE_A ] ){   }
	else if ( keys[ SDL_SCANCODE_D ] ){   }

    if      ( keys[ SDL_SCANCODE_W ] ){  }
	else if ( keys[ SDL_SCANCODE_S ] ){   }

    if      ( keys[ SDL_SCANCODE_E ] ){ }
	else if ( keys[ SDL_SCANCODE_Q ] ){}

    //if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[SDL_SCANCODE_D]||keys[SDL_SCANCODE_E]||keys[SDL_SCANCODE_Q] ){
    if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_E]||keys[SDL_SCANCODE_Q] ){
        controler.bActive = false;
        Mat3d rot;
        //rot.setT(myCraft->rotMat);
        rot.set(myCraft->rotMat);
        controler.goalDir = rot.c;
        controler.goalUp  = rot.b;
    }else{
        controler.bActive = true;
    };


    Mat3d rot; rot.setT(myCraft->rotMat);
	//world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
}
*/


void MusicVisualizerGUI::eventHandling( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                   // mouseSteer = true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                  //  mouseSteer = false;
                 //   pilot->resetSteer();
                    break;
            }
            break;
    }

};

MusicVisualizerGUI * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new MusicVisualizerGUI( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}

