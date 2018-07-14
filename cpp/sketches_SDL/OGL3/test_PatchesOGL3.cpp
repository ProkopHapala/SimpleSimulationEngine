

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

#include "Solids.h"

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
#include "GL3Utils.h"

#include "IO_utils.h"


class HorizontTestApp : public AppSDL2OGL3, public SceneOGL3 { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    float camDist = 100.0;
    float screenDist = 100.0;

    Shader *sh1=0,*shTx=0;
    GLMesh *oglQuad1=0, *oglTri1=0, *oglHalfHex1=0;

	//virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );

    virtual void draw( Camera& cam );

	HorizontTestApp(int W, int H);

};

HorizontTestApp::HorizontTestApp(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){

    //DEBUG_mesh = new GLMeshBuilder();
    //DEBUG_VIEW_INIT()

    for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    //shDebug=new Shader();
    //shDebug->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    //shDebug->getDefaultUniformLocation();

    sh1=new Shader();
    sh1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    sh1->getDefaultUniformLocation();

    shTx=new Shader();
    shTx->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTx->getDefaultUniformLocation();

    sh1->use();
    glUniform3f(sh1->getUloc("lightColor"   ), 0.5,0.45,0.4 );
    glUniform3f(sh1->getUloc("diffuseColor" ), 1.0,1.0,1.0  );
    glUniform3f(sh1->getUloc("ambientColor" ), 0.2,0.25,0.3 );
    glUniform3f(sh1->getUloc("specularColor"), 2.0,2.0,2.0  );
    glUniform3f(sh1->getUloc("lightPos"     ), 10.0,-10.0,-10.0 );
    //glUniform3f(sh1->getUloc("lightPos"     ), 10.0,10.0,10.0 );

    oglQuad1    = glQuadGrid    ( {10,10}, true );  //oglQuad1->draw_mode = GL_LINE_STRIP;
    oglTri1     = glTriangleGrid( 10    , true  );  //oglTri1 ->draw_mode = GL_LINE_STRIP;
    oglHalfHex1 = glHalfHexGrid( {10,15}, true );   //oglHalfHex1 ->draw_mode = GL_LINE_STRIP;

    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0;
    cam.zmax = 1000.0;
    cam.zoom = 5.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;


};

//void HorizontTestApp::update(){
//    AppSDL2OGL3::update();
//}

void HorizontTestApp::draw( Camera& cam ){

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0,0,screen->WIDTH,screen->HEIGHT);
    glScissor (0,0,screen->WIDTH,screen->HEIGHT);

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    sh1->use();
    //cam.lookAt( Vec3dZ, 20.0 );
    setCamera( *sh1, cam );
    //setCameraOrtho(*sh1,cam);
    //sh1->setModelPose( myCraft->pos, myCraft->rotMat );
    sh1->setModelPoseT( Vec3dZero, Mat3dIdentity );

    GLuint ucolor = sh1->getUloc("baseColor");
    glUniform4f( ucolor, 0.0, 0.0, 0.0, 1.0 );
    //oglQuad1->draw_mode = GL_LINES;

    //oglQuad1->draw();
    //oglTri1->draw();
    oglHalfHex1->draw();

};

void HorizontTestApp::eventHandling( const SDL_Event& event  ){
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

HorizontTestApp * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new HorizontTestApp( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}



